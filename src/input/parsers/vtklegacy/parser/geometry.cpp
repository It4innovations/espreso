
#include "geometry.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/parser.h"
#include "config/ecf/input/input.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "mesh/element.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"
#include "input/parsers/asciiparser.h"
#include "input/parsers/mixedelementsparser.h"

#include <numeric>

using namespace espreso;

static int _dimension(int vtktype)
{
	switch (vtktype) {
	case  1: return 0;
	case  3: return 1;
	case 21: return 1;
	case  9: return 2;
	case 23: return 2;
	case  5: return 2;
	case 22: return 2;
	case 10: return 3;
	case 24: return 3;
	case 14: return 3;
	case 27: return 3;
	case 13: return 3;
	case 26: return 3;
	case 12: return 3;
	case 25: return 3;
	}
	return -1;
}

static int _etype(int dimension, int esize, const std::string &file)
{
	switch (dimension) {
	case 0:
		break; // points are not included into elements
	case 1:
		switch (esize) {
		case 2: return (int)Element::CODE::LINE2;
		case 3: return (int)Element::CODE::LINE3;
		} break;
	case 2:
		switch (esize) {
		case 3: return (int)Element::CODE::TRIANGLE3;
		case 4: return (int)Element::CODE::SQUARE4;
		case 6: return (int)Element::CODE::TRIANGLE6;
		case 8: return (int)Element::CODE::SQUARE8;
		} break;
	case 3:
		switch (esize) {
		case  4: return (int)Element::CODE::TETRA4;
		case  5: return (int)Element::CODE::PYRAMID5;
		case  6: return (int)Element::CODE::PRISMA6;
		case  8: return (int)Element::CODE::HEXA8;
		case 10: return (int)Element::CODE::TETRA10;
		case 13: return (int)Element::CODE::PYRAMID13;
		case 15: return (int)Element::CODE::PRISMA15;
		case 20: return (int)Element::CODE::HEXA20;
		} break;
	}
	eslog::error("VTK Legacy parser: unrecognized element type (dim=%d,esize=%d) in file='%s'.\n", dimension, esize, file.c_str());
	return -1;
}

VTKLegacyGeometry::Keyword::Keyword(InputFilePack &pack, const char *c)
: begin(0), end(0)
{
	fileindex = pack.fileindex;
	rank = info::mpi::rank;
	offset = pack.distribution[info::mpi::rank] + (c - pack.begin);
}

VTKLegacyGeometry::Header::Header(InputFilePack &pack, const char *c)
: Keyword(pack, c), format(Format::UNKNOWN), dataset(DataSet::UNKNOWN)
{
	while (*c++ != '\n'); // version
	while (*c++ != '\n'); // header
	if (StringCompare::caseInsensitiveEq(std::string(c, c + 5), "ASCII")) {
		format = Format::ASCII;
	}
	if (StringCompare::caseInsensitiveEq(std::string(c, c + 5), "BINARY")) {
		format = Format::BINARY;
	}
	while (*c++ != '\n'); // format
	if (StringCompare::caseInsensitiveEq(std::string(c, c + 7), "DATASET")) {
		if (StringCompare::caseInsensitiveEq(std::string(c + 8, c + 8 + 17), "UNSTRUCTURED_GRID")) {
			dataset = DataSet::UNSTRUCTURED_GRID;
		}
	}
}

VTKLegacyGeometry::Points::Points(InputFilePack &pack, const char *c)
: Keyword(pack, c)
{
	const char *_c = c;
	nn = strtol(c + 7, NULL, 0);
	c += 7; // POINTS
	c += ASCIIParser::keyend(c); // nn
	c += ASCIIParser::keyend(c); // format
	begin = offset + (c - _c);
}

VTKLegacyGeometry::Cells::Cells(InputFilePack &pack, const char *c)
: Keyword(pack, c)
{
	char *next;
	const char *_c = c;
	ne = strtol(c + 6, &next, 0);
	c = next;
	size = strtol(c, &next, 0);
	begin = offset + (next - _c);
}

VTKLegacyGeometry::CellTypes::CellTypes(InputFilePack &pack, const char *c)
: Keyword(pack, c)
{
	char *next;
	ne = strtol(c + 11, &next, 0);
	begin = offset + (next - c);
}

VTKLegacyGeometry::Data::Data(InputFilePack &pack, DataSource source, const char *c)
: Keyword(pack, c), source(source)
{
	// TODO: finish implementation
	begin = offset + (source == DataSource::CELLS ? 9 : 10);
}

VTKLegacyGeometry::VTKLegacyGeometry(InputFilePack &pack)
: _pack(pack)
{

}

void VTKLegacyGeometry::scan()
{
	header();
	scanASCII();
	scanBinary();
}

void VTKLegacyGeometry::parse(MeshBuilder &mesh, const std::vector<std::string> &names)
{
	parseASCII(mesh, names);
	parseBinary(mesh, names);
}

void VTKLegacyGeometry::header()
{
	while (_pack.next()) {
		if (_pack.distribution[info::mpi::rank] == 0 && _pack.distribution[info::mpi::rank + 1] != 0) {
			_header.push_back({ _pack, _pack.begin });
			if (_header.back().format == Format::UNKNOWN) {
				eslog::error("VTK Legacy parser: file '%s' has unknown VTK file format.\n", _pack.files[_pack.fileindex]->name.c_str());
			}
			if (_header.back().dataset == DataSet::UNKNOWN) {
				eslog::error("VTK Legacy parser: file '%s' unsupported DATASET TYPE.\n", _pack.files[_pack.fileindex]->name.c_str());
			}
		}
	}
	DistributedScanner parser;
	parser.synchronize(_header);
}

template <typename TKeyword>
static void addoffset(std::vector<std::vector<size_t> > &offsets, std::vector<TKeyword> &keywords)
{
	for (auto keyword = keywords.begin(); keyword != keywords.end(); ++keyword) {
		offsets[keyword->fileindex].push_back(keyword->offset);
	}
}

template <typename TKeyword, typename ...TOther>
static void addoffset(std::vector<std::vector<size_t> > &offsets, std::vector<TKeyword> &keywords, TOther& ...other)
{
	addoffset(offsets, keywords);
	addoffset(offsets, other...);
}

template <typename TKeyword>
static void setend(std::vector<std::vector<size_t> > &offsets, std::vector<TKeyword> &keywords)
{
	for (auto keyword = keywords.begin(); keyword != keywords.end(); ++keyword) {
		keyword->end = *std::lower_bound(offsets[keyword->fileindex].begin(), offsets[keyword->fileindex].end(), keyword->begin);
	}
}

template <typename TKeyword, typename ...TOther>
static void setend(std::vector<std::vector<size_t> > &offsets, std::vector<TKeyword> &keywords, TOther& ...other)
{
	setend(offsets, keywords);
	setend(offsets, other...);
}

void VTKLegacyGeometry::scanBinary()
{
//	DistributedParser parser;

	while (_pack.next()) {
		if (_header[_pack.fileindex].format == Format::BINARY) {
			eslog::error("VTK Legacy parser: not implemented scanning of binary format.\n");
			// parser.scan(_pack);
		}
	}
	// parser.synchronize(_points, _cells, _cellTypes, _pointData, _cellData);
}

void VTKLegacyGeometry::parseBinary(MeshBuilder &mesh, const std::vector<std::string> &names)
{
//	while (_pack.next()) {
//		if (_header[_pack.fileindex].format == Format::BINARY) {
//			// parse
//		}
//	}
}

void VTKLegacyGeometry::scanASCII()
{
	// TODO: fully case insensitive ???, e.g. Points, POints, etc..
	DistributedScanner scanner;

	scanner.add({ "points", "POINTS" }, [&] (const char *c) { _points.push_back({_pack, c}); }, [&] (const char *c) { return 2 * Points(_pack, c).nn; });
	scanner.add({ "cells", "CELLS" }, [&] (const char *c) { _cells.push_back({_pack, c}); }, [&] (const char *c) { return 2 * Cells(_pack, c).size; });
	scanner.add({ "cell_types", "CELL_TYPES" }, [&] (const char *c) { _cellTypes.push_back({_pack, c}); }, [&] (const char *c) { return 2 * CellTypes(_pack, c).ne; });
	scanner.add({ "point_data", "POINT_DATA" }, [&] (const char *c) { _pointData.push_back({_pack, DataSource::POINTS, c}); });
	scanner.add({ "cell_data", "CELL_DATA" }, [&] (const char *c) { _cellData.push_back({_pack, DataSource::CELLS, c}); });

	while (_pack.next()) {
		if (_header[_pack.fileindex].format == Format::ASCII) {
			scanner.align(_pack, " \n");
			scanner.scanlines(_pack);
		}
	}

	scanner.synchronize(_points, _cells, _cellTypes, _pointData, _cellData);

	std::sort(_points.begin(), _points.end(), [] (Points &p1, Points &p2) { return p1.fileindex < p2.fileindex; });
	std::sort(_cells.begin(), _cells.end(), [] (Cells &c1, Cells &c2) { return c1.fileindex < c2.fileindex; });
	std::sort(_cellTypes.begin(), _cellTypes.end(), [] (CellTypes &t1, CellTypes &t2) { return t1.fileindex < t2.fileindex; });

	std::vector<std::vector<size_t> > offsets(_pack.files.size());
	addoffset(offsets, _points, _cells, _cellTypes, _pointData, _cellData);
	for (size_t i = 0; i < offsets.size(); ++i) {
		offsets[i].push_back(_pack.files[i]->distribution.back());
		std::sort(offsets[i].begin(), offsets[i].end());
	}
	setend(offsets, _points, _cells, _cellTypes, _pointData, _cellData);
}

void VTKLegacyGeometry::parseASCII(MeshBuilder &mesh, const std::vector<std::string> &names)
{
	std::vector<std::vector<double> > points(_pack.files.size());
	std::vector<std::vector<esint> > cells(_pack.files.size()), celltypes(_pack.files.size());

	std::vector<size_t> npoints(_pack.files.size());
	std::vector<int> mindim(_pack.files.size(), 4), maxdim(_pack.files.size()), _mindim(_pack.files.size(), 4), _maxdim(_pack.files.size());

	while (_pack.next()) {
		// only mandatory parts
		// data are currently skipped
		ASCIIParser::parse(points[_pack.fileindex], _pack, _points[_pack.fileindex].begin, _points[_pack.fileindex].end);
		ASCIIParser::parse(cells[_pack.fileindex], _pack, _cells[_pack.fileindex].begin, _cells[_pack.fileindex].end);
		ASCIIParser::parse(celltypes[_pack.fileindex], _pack, _cellTypes[_pack.fileindex].begin, _cellTypes[_pack.fileindex].end);

		npoints[_pack.fileindex] = points[_pack.fileindex].size();
		for (size_t i = 0; i < celltypes[_pack.fileindex].size(); ++i) {
			int d = _dimension(celltypes[_pack.fileindex][i]);
			_mindim[_pack.fileindex] = std::min(_mindim[_pack.fileindex], d);
			_maxdim[_pack.fileindex] = std::max(_maxdim[_pack.fileindex], d);
		}
	}

	std::vector<size_t> dummy(_pack.files.size());
	Communication::exscan(dummy, npoints);
	Communication::allReduce(_mindim.data(), mindim.data(), _pack.files.size(), MPI_INT, MPI_MIN);
	Communication::allReduce(_maxdim.data(), maxdim.data(), _pack.files.size(), MPI_INT, MPI_MAX);

	while (_pack.next()) {
		if (mindim[_pack.fileindex] != maxdim[_pack.fileindex]) {
			eslog::globalerror("VTK Legacy parser: not implemented parsing of a file with various elements dimension: '%s'.\n", _pack.files[_pack.fileindex]->name.c_str());
		}
	}

	MixedElementsParser mixedparser;
	for (size_t i = 0; i < cells.size(); ++i) {
		mixedparser.add(cells[i].data(), cells[i].size());
	}
	mixedparser.parse([&] (size_t index, esint id) {
		switch (mindim[index]) {
		case 0:
			if (id == 1) { return id; }
			break;
		case 1:
			if (id == 2 || id == 3) { return id; }
			break;
		case 2:
			if (id == 3 || id == 4 || id == 6 || id == 8) { return id; }
			break;
		case 3:
			if (id == 4 || id == 5 || id == 6 || id == 8 || id == 10 || id == 13 || id == 15 || id == 20) { return id; }
			break;
		}
		return esint{};
	});

	for (size_t i = 0; i < mixedparser.invalid.size(); ++i) {
		if (mixedparser.invalid[i] != info::mpi::size) {
			// it happens in really rare cases
			// with the increasing topology size, the probability decreases
			eslog::warning("VTK Legacy parser: synchronization of region ''.\n", names[i].c_str());
		}
	}

	size_t csum = 0, esum = 0;
	while (_pack.next()) {
		size_t pbegin = (3 - (npoints[_pack.fileindex] % 3)) % 3;
		size_t pmissing = (3 - ((points[_pack.fileindex].size() - pbegin) % 3)) % 3;
		ASCIIParser::addmore(points[_pack.fileindex], _pack, pmissing, _points[_pack.fileindex].end);
		ASCIIParser::addmore(cells[_pack.fileindex], _pack, mixedparser.missing[_pack.fileindex], _cells[_pack.fileindex].end);

		csum += (points[_pack.fileindex].size() - pbegin) / 3;
		esum += cells[_pack.fileindex].size();
	}

	mesh.nIDs.reserve(csum);
	mesh.coordinates.reserve(csum);
	mesh.esize.reserve(esum / 2);
	mesh.etype.reserve(esum / 2);
	mesh.eIDs.reserve(esum / 2);
	mesh.enodes.reserve(esum / 2);

	size_t nidoffset = 0, eidoffset = 0;
	while (_pack.next()) {
		{ // fill coordinates
			size_t pbegin = (3 - (npoints[_pack.fileindex] % 3)) % 3;
			size_t csize = (points[_pack.fileindex].size() - pbegin) / 3;
			size_t noffset = (npoints[_pack.fileindex] + pbegin) / 3;
			for (size_t i = 0; i < csize; ++i) {
				mesh.nIDs.push_back(nidoffset + i + noffset);
				mesh.coordinates.push_back(Point(
						points[_pack.fileindex][3 * i + pbegin + 0],
						points[_pack.fileindex][3 * i + pbegin + 1],
						points[_pack.fileindex][3 * i + pbegin + 2]));
			}
		}
		{ // fill elements ids, size, and nodes
			std::vector<esint> ids;
			int nelements = 0;
			for (size_t i = mixedparser.first[_pack.fileindex]; i < cells[_pack.fileindex].size(); i += cells[_pack.fileindex][i] + 1) {
				if (cells[_pack.fileindex][i] == 1) {
					ids.push_back(cells[_pack.fileindex][i + 1] + nidoffset);
				} else {
					mesh.esize.push_back(cells[_pack.fileindex][i]);
					mesh.etype.push_back(_etype(maxdim[_pack.fileindex], mesh.esize.back(), names[_pack.fileindex]));
					for (int n = 0; n < mesh.esize.back(); ++n) {
						mesh.enodes.push_back(cells[_pack.fileindex][i + 1 + n] + nidoffset);
					}
					ids.push_back(eidoffset + mixedparser.offset[_pack.fileindex] + nelements++);
				}
			}

			if (maxdim[_pack.fileindex]) {
				mesh.eregions[names[_pack.fileindex]] = ids;
				mesh.eIDs.insert(mesh.eIDs.end(), ids.begin(), ids.end());
			} else {
				mesh.nregions[names[_pack.fileindex]] = ids;
			}
		}
		nidoffset += _points[_pack.fileindex].nn;
		eidoffset += _cells[_pack.fileindex].ne;
	}
}

