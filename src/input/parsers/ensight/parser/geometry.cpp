
#include "geometry.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "mesh/element.h"
#include "input/input.h"
#include "input/parsers/fileblock.h"
#include "input/parsers/distributedscanner.h"

#include <numeric>

using namespace espreso;

EnsightGeometry::EnsightGeometry(InputFilePack &geofile, OrderedMeshDatabase &database)
: _geofile(geofile), _database(database)
{

}

void EnsightGeometry::scan()
{
	// root always has header
	if (info::mpi::rank == 0) {
		const char *ndesc = _geofile.begin, *edesc = _geofile.begin;
		if (memcmp(_geofile.begin, "C Binary", 8) == 0) {
			_keywords.header.format = EnsightKeywords::Format::BINARY;
			ndesc += 240, edesc += 320;
		} else {
			_keywords.header.format = EnsightKeywords::Format::ASCII;
			while (*ndesc++ != '\n');
			while (*ndesc++ != '\n');
			edesc = ndesc;
			while (*edesc++ != '\n');
		}
		if (DistributedScanner::check(ndesc, "node id ")) {
			if (DistributedScanner::check(ndesc + 8, "off"))    { _keywords.header.nodeIDs = EnsightKeywords::IDs::OFF; }
			if (DistributedScanner::check(ndesc + 8, "given"))  { _keywords.header.nodeIDs = EnsightKeywords::IDs::GIVEN; }
			if (DistributedScanner::check(ndesc + 8, "assign")) { _keywords.header.nodeIDs = EnsightKeywords::IDs::ASSIGN; }
			if (DistributedScanner::check(ndesc + 8, "ignore")) { _keywords.header.nodeIDs = EnsightKeywords::IDs::INGNORE; }
		}
		if (DistributedScanner::check(edesc, "element id ")) {
			if (DistributedScanner::check(edesc + 11, "off"))    { _keywords.header.elementIDs = EnsightKeywords::IDs::OFF; }
			if (DistributedScanner::check(edesc + 11, "given"))  { _keywords.header.elementIDs = EnsightKeywords::IDs::GIVEN; }
			if (DistributedScanner::check(edesc + 11, "assign")) { _keywords.header.elementIDs = EnsightKeywords::IDs::ASSIGN; }
			if (DistributedScanner::check(edesc + 11, "ignore")) { _keywords.header.elementIDs = EnsightKeywords::IDs::INGNORE; }
		}
	}

	Communication::broadcast(&_keywords.header, sizeof(EnsightKeywords::Header), MPI_BYTE, 0);

	if (_keywords.header.format == EnsightKeywords::Format::UNKNOWN || _keywords.header.nodeIDs == EnsightKeywords::IDs::UNKNOWN || _keywords.header.elementIDs == EnsightKeywords::IDs::UNKNOWN) {
		eslog::globalerror("MESIO internal error: cannot recognize EnSight Gold geometry format.\n");
	}

	DistributedScanner scanner;
	if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
		fillScanner<EnsightASCIIGeometryKeywordParser>(_geofile, scanner, _keywords);
	} else {
		fillScanner<EnsightBinaryGeometryKeywordParser>(_geofile, scanner, _keywords);
	}

	scanner.scan(_geofile);
	scanner.synchronize(_keywords.parts, _keywords.coordinates, _keywords.elements);
}

void EnsightGeometry::parse()
{
	size_t csize, esize, esuffix;
	if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
		csize = 13;
		esize = 10;
		esuffix = 1;
	} else {
		csize = sizeof(float);
		esize = sizeof(int);
		esuffix = 0;
	}

	// skip ids
	if (_keywords.header.nodeIDs == EnsightKeywords::IDs::GIVEN || _keywords.header.nodeIDs == EnsightKeywords::IDs::INGNORE) {
		for (size_t i = 0; i < _keywords.coordinates.size(); ++i) {
			_keywords.coordinates[i].offset += _keywords.coordinates[i].nn * (esize + esuffix);
		}
	}
	if (_keywords.header.elementIDs == EnsightKeywords::IDs::GIVEN || _keywords.header.elementIDs == EnsightKeywords::IDs::INGNORE) {
		for (size_t i = 0; i < _keywords.elements.size(); ++i) {
			_keywords.elements[i].offset += _keywords.elements[i].ne * (esize + esuffix);
		}
	}

	esint coffset = 0, eoffset = 0;
	for (size_t p = 0, e = 0; p < _keywords.parts.size(); ++p) {
		std::string name = _keywords.parts[p].getName();

		{ // coordinate block
			std::vector<float> coordinates;
			FileBlock block(_geofile, _keywords.coordinates[p].offset, 3 * csize * _keywords.coordinates[p].nn, csize, info::mpi::rank);
			if (block.size) {
				if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
					coordinates.reserve(block.size / csize);
					for (const char *cc = _geofile.begin + block.begin; cc < _geofile.begin + block.end; cc += csize) {
						coordinates.push_back(atof(cc));
					}
				} else {
					coordinates.resize(block.size / csize);
					memcpy(coordinates.data(), _geofile.begin + block.begin, block.size);
				}

				int start = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _keywords.coordinates[p].offset) - _geofile.distribution.begin() - 1;
				int end = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _keywords.coordinates[p].offset + 3 * csize * _keywords.coordinates[p].nn) - _geofile.distribution.begin() - 1;

				if (start == info::mpi::rank) {
					std::vector<float, initless_allocator<float> > rbuffer(3 * _keywords.coordinates[p].nn);
					std::copy(coordinates.begin(), coordinates.end(), rbuffer.begin());
					std::vector<MPI_Request> req(end - start);
					for (int r = start + 1; r <= end; ++r) {
						FileBlock rblock(_geofile, _keywords.coordinates[p].offset, csize * _keywords.coordinates[p].nn, csize, r);
						MPI_Irecv(rbuffer.data() + rblock.prevsize / csize, rblock.size / csize, MPI_FLOAT, r, 0, info::mpi::comm, req.data() + (r - start - 1));
					}
					_database.noffset.reserve(_database.coordinates.size() + _keywords.coordinates[p].nn);
					_database.coordinates.reserve(_database.coordinates.size() + _keywords.coordinates[p].nn);
					MPI_Waitall(end - start, req.data(), MPI_STATUSES_IGNORE);
					for (int n = 0; n < _keywords.coordinates[p].nn; ++n) {
						_database.noffset.push_back(coffset + n);
						_database.coordinates.push_back(_Point<esfloat>(rbuffer[n + _keywords.coordinates[p].nn * 0], rbuffer[n + _keywords.coordinates[p].nn * 1], rbuffer[n + _keywords.coordinates[p].nn * 2]));
					}
				} else {
					MPI_Send(coordinates.data(), coordinates.size(), MPI_FLOAT, start, 0, info::mpi::comm);
				}
			}
		}

		for (; e < _keywords.elements.size() && (p + 1 == _keywords.parts.size() || _keywords.elements[e].offset < _keywords.coordinates[p + 1].offset); ++e) {
			std::vector<int> elements;
			int enodes = _keywords.elements[e].getSize();
			FileBlock block(_geofile, _keywords.elements[e].offset, (enodes * esize + esuffix) * _keywords.elements[e].ne, enodes * esize + esuffix, info::mpi::rank);
			if (block.size) {
				if (_keywords.header.format == EnsightKeywords::Format::ASCII) {
					elements.reserve(enodes * ((block.end - block.begin) / (enodes * esize + esuffix)));
					for (const char *d = _geofile.begin + block.begin; d < _geofile.begin + block.end; d += enodes * esize + esuffix) {
						for (int n = 0; n < enodes; ++n) {
							elements.push_back(atol(d + 10 * n));
						}
					}
				} else {
					elements.resize(enodes * block.size / (enodes * esize + esuffix));
					memcpy(elements.data(), _geofile.begin + block.begin, block.size);
				}

				if (_keywords.elements[e].getCode() != Element::CODE::POINT1) {
					esint boffset = block.prevsize / (enodes * esize + esuffix);
					esint bsize = block.size / (enodes * esize + esuffix);
					_database.etype.resize(_database.etype.size() + bsize, (int)_keywords.elements[e].getCode());
					_database.esize.resize(_database.esize.size() + bsize, enodes);
					_database.enodes.reserve(_database.enodes.size() + elements.size());
					for (esint i = 0; i < bsize; ++i) {
						_database.eoffset.push_back(eoffset + boffset + i);
					}
					for (size_t n = 0; n < elements.size(); ++n) {
						_database.enodes.push_back(elements[n] + coffset);
					}
				}
			}

			if (_keywords.elements[e].getCode() == Element::CODE::POINT1) {
				_database.nregions.push_back(OrderedMeshDatabase::Region{ name, coffset, coffset + _keywords.coordinates[p].nn });
			} else {
				if (_database.eregions.empty() || _database.eregions.back().name.compare(name)) {
					_database.eregions.push_back(OrderedMeshDatabase::Region{ name, eoffset, eoffset + _keywords.elements[e].ne });
				} else {
					_database.eregions.back().end += _keywords.elements[e].ne;
				}
				eoffset += _keywords.elements[e].ne;
			}
		}
		coffset += _keywords.coordinates[p].nn;
	}
}
