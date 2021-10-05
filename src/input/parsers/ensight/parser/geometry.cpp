
#include "geometry.h"
#include "wrappers/mpi/communication.h"
#include "esinfo/eslog.h"
#include "mesh/element.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"

#include <numeric>

using namespace espreso;
//test
EnsightGeometry::EnsightGeometry(InputFilePack &geofile)
: _geofile(geofile), _header{ Format::UNKNOWN, IDs::UNKNOWN, IDs::UNKNOWN }
{

}

void EnsightGeometry::scan()
{
	header();

	if (_header.format == Format::ASCII) {
		scanASCII();
	} else {
		scanBinary();
	}
}

void EnsightGeometry::parse(MeshBuilder &mesh)
{
	if (_header.format == Format::ASCII) {
		parseASCII(mesh);
	} else {
		parseBinary(mesh);
	}
}

void EnsightGeometry::header()
{
	// root always has header
	if (info::mpi::rank == 0) {
		if (memcmp(_geofile.begin, "C Binary", 8) == 0) {
			_header.format = Format::BINARY;

			size_t ndesc = 240, edesc = 320;
			if (DistributedScanner::check(_geofile.begin + ndesc, "node id ")) {
				if (DistributedScanner::check(_geofile.begin + ndesc + 8, "off")) { _header.nodeIDs = IDs::OFF; }
				if (DistributedScanner::check(_geofile.begin + ndesc + 8, "given")) { _header.nodeIDs = IDs::GIVEN; }
				if (DistributedScanner::check(_geofile.begin + ndesc + 8, "assign")) { _header.nodeIDs = IDs::ASSIGN; }
				if (DistributedScanner::check(_geofile.begin + ndesc + 8, "ignore")) { _header.nodeIDs = IDs::INGNORE; }
			}
			if (DistributedScanner::check(_geofile.begin + edesc, "element id ")) {
				if (DistributedScanner::check(_geofile.begin + edesc + 11, "off")) { _header.elementIDs = IDs::OFF; }
				if (DistributedScanner::check(_geofile.begin + edesc + 11, "given")) { _header.elementIDs = IDs::GIVEN; }
				if (DistributedScanner::check(_geofile.begin + edesc + 11, "assign")) { _header.elementIDs = IDs::ASSIGN; }
				if (DistributedScanner::check(_geofile.begin + edesc + 11, "ignore")) { _header.elementIDs = IDs::INGNORE; }
			}
		} else {
			_header.format = Format::ASCII;
		}
	}

	Communication::broadcast(&_header, sizeof(Header), MPI_BYTE, 0);

	/*if (_header.format == Format::UNKNOWN || _header.nodeIDs == IDs::UNKNOWN || _header.elementIDs == IDs::UNKNOWN) {
		eslog::globalerror("MESIO internal error: cannot recognize EnSight Gold geometry format.\n");
	}*/
}

void EnsightGeometry::scanBinary()
{
	DistributedScanner parser;

	auto getsize = [] (const void *data) {
		int size;
		memcpy(&size, data, sizeof(int));
		return size;
	};

	auto addcoordinates = [&] (const char *c) {
		_coordinates.push_back(Coordinates(_geofile.distribution[info::mpi::rank] + (c - _geofile.begin) + 80 + sizeof(int), getsize(c + 80)));
	};
	auto skipcoordinates = [&] (const char *c) -> size_t {
		int nn; memcpy(&nn, c + 80, sizeof(int));
		return 3 * nn * sizeof(float);
	};

	auto addelement = [&] (const char *c, Elements::Type type) {
		size_t desc = 80 + sizeof(int), offset = _geofile.distribution[info::mpi::rank] + (c - _geofile.begin);
		int nn; memcpy(&nn, c + 80, sizeof(int));
		_elements.push_back(Elements(type, offset + desc, nn));
	};

	auto skipelements = [&] (const char *c, int enodes) -> size_t {
		int nn; memcpy(&nn, c + 80, sizeof(int));
		return enodes * nn * sizeof(int);
	};

	parser.add("part"       , [&] (const char *c) { _parts.insert(_parts.end(), c + 80 + sizeof(int), c + 80 + sizeof(int) + 80); });
	parser.add("coordinates", [&] (const char *c) { addcoordinates(c); }                       , [&] (const char *c) { return skipcoordinates(c); });

	parser.add("point"      , [&] (const char *c) { addelement(c, Elements::Type::POINT); }    , [&] (const char *c) { return skipelements(c,  1); });
	parser.add("bar2"       , [&] (const char *c) { addelement(c, Elements::Type::BAR2); }     , [&] (const char *c) { return skipelements(c,  2); });
	parser.add("bar3"       , [&] (const char *c) { addelement(c, Elements::Type::BAR3); }     , [&] (const char *c) { return skipelements(c,  3); });
	parser.add("tria3"      , [&] (const char *c) { addelement(c, Elements::Type::TRIA3); }    , [&] (const char *c) { return skipelements(c,  3); });
	parser.add("tria6"      , [&] (const char *c) { addelement(c, Elements::Type::TRIA6); }    , [&] (const char *c) { return skipelements(c,  6); });
	parser.add("quad4"      , [&] (const char *c) { addelement(c, Elements::Type::QUAD4); }    , [&] (const char *c) { return skipelements(c,  4); });
	parser.add("quad8"      , [&] (const char *c) { addelement(c, Elements::Type::QUAD8); }    , [&] (const char *c) { return skipelements(c,  8); });
	parser.add("tetra4"     , [&] (const char *c) { addelement(c, Elements::Type::TETRA4); }   , [&] (const char *c) { return skipelements(c,  4); });
	parser.add("tetra10"    , [&] (const char *c) { addelement(c, Elements::Type::TETRA10); }  , [&] (const char *c) { return skipelements(c, 10); });
	parser.add("pyramid5"   , [&] (const char *c) { addelement(c, Elements::Type::PYRAMID5); } , [&] (const char *c) { return skipelements(c,  5); });
	parser.add("pyramid13"  , [&] (const char *c) { addelement(c, Elements::Type::PYRAMID13); }, [&] (const char *c) { return skipelements(c, 13); });
	parser.add("penta6"     , [&] (const char *c) { addelement(c, Elements::Type::PENTA6); }   , [&] (const char *c) { return skipelements(c,  6); });
	parser.add("penta15"    , [&] (const char *c) { addelement(c, Elements::Type::PENTA15); }  , [&] (const char *c) { return skipelements(c, 15); });
	parser.add("hexa8"      , [&] (const char *c) { addelement(c, Elements::Type::HEXA8); }    , [&] (const char *c) { return skipelements(c,  8); });
	parser.add("hexa20"     , [&] (const char *c) { addelement(c, Elements::Type::HEXA20); }   , [&] (const char *c) { return skipelements(c, 20); });

	parser.add("nsided"     , [&] (const char *c) { eslog::error("Ensight parser error: ESPRESO does not support nsided elements.\n"); });
	parser.add("nfaced"     , [&] (const char *c) { eslog::error("Ensight parser error: ESPRESO does not support nfaced elements.\n"); });

	parser.scan(_geofile);
	parser.synchronize(_parts, _coordinates, _elements);
}

void EnsightGeometry::parseBinary(MeshBuilder &mesh)
{
	// skip ids: currently we always generate them
	if (_header.nodeIDs == IDs::GIVEN || _header.nodeIDs == IDs::INGNORE) {
		for (size_t i = 0; i < _coordinates.size(); ++i) {
			_coordinates[i].offset += _coordinates[i].nn * sizeof(int);
		}
	}
	if (_header.elementIDs == IDs::GIVEN || _header.elementIDs == IDs::INGNORE) {
		for (size_t i = 0; i < _elements.size(); ++i) {
			_elements[i].offset += _elements[i].ne * sizeof(int);
		}
	}

	auto round = [] (size_t &current, size_t &start, size_t chunk) {
		if ((current - start) % chunk != 0) {
			current += chunk - (current - start) % chunk;
		}
	};

	struct eoffsetinfo {
		std::string name;
		int offset;
		int size;
		int totalsize;
	};

	std::vector<eoffsetinfo> einfo;
	size_t parts = _parts.size() / 80;
	esint cidoffset = 0;

	for (size_t p = 0, i = 0; p < parts; ++p) {
		auto begin = _parts.begin() + p * 80;
		auto end = _parts.begin() + p * 80;
		while (*end != 0) { ++end; }
		std::string name = std::string(begin, end);

		int firstrank = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[p].offset) - _geofile.distribution.begin() - 1;
		int lastrank = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[p].offset + 3 * sizeof(float) * _coordinates[p].nn) - _geofile.distribution.begin() - 1;

		size_t cbegin = std::max(_coordinates[p].offset, _geofile.distribution[info::mpi::rank]);
		size_t cend = std::min(_coordinates[p].offset + 3 * sizeof(float) * _coordinates[p].nn, _geofile.distribution[info::mpi::rank + 1]);
		if (cbegin < cend) {
			// read data: x x x x x...; y y y y y y...; z z z z z z...
			round(cbegin, _coordinates[p].offset, sizeof(float));
			round(cend, _coordinates[p].offset, sizeof(float));
			std::vector<float> coordinates((cend - cbegin) / sizeof(float));
			memcpy(coordinates.data(), _geofile.begin + cbegin - _geofile.distribution[info::mpi::rank], cend - cbegin);

			std::vector<size_t> cdistribution = tarray<size_t>::distribute(lastrank - firstrank + 1, _coordinates[p].nn);
			if (firstrank != lastrank) { // group data: x y z; x y z; x y z; ...
				cbegin = (cbegin - _coordinates[p].offset) / sizeof(float); // from char to float offset
				cend   = (cend   - _coordinates[p].offset) / sizeof(float); // from char to float offset
				std::vector<int> sBuffer, rBuffer;
				for (int r = firstrank; r <= lastrank; ++r) {
					size_t prevsize = sBuffer.size();
					sBuffer.push_back(0); // total size
					sBuffer.push_back(r); // target
					sBuffer.push_back(info::mpi::rank); // source
					sBuffer.push_back(0); // x size
					sBuffer.push_back(0); // y size
					sBuffer.push_back(0); // z size

					for (int d = 0; d < 3; ++d) {
						size_t rbegin = std::max(d * _coordinates[p].nn + cdistribution[r - firstrank], cbegin);
						size_t rend   = std::min(d * _coordinates[p].nn + cdistribution[r - firstrank + 1], cend);
						if (rend < rbegin) {
							rend = rbegin;
						}
						sBuffer[prevsize + 3 + d] = rend - rbegin;
						if (rbegin != rend) {
							sBuffer.insert(sBuffer.end(), reinterpret_cast<const int*>(coordinates.data() + rbegin - cbegin), reinterpret_cast<const int*>(coordinates.data() + rend - cbegin));
						}
					}
					sBuffer[prevsize] = sBuffer.size() - prevsize;
				}

				if (!Communication::allToAllWithDataSizeAndTarget(sBuffer, rBuffer, firstrank, lastrank + 1)) {
					eslog::error("Ensight Gold parser: cannot exchange coordinate data.\n");
				}

				coordinates.resize(3 * (cdistribution[info::mpi::rank - firstrank + 1] - cdistribution[info::mpi::rank - firstrank]));
				std::vector<int> alloffsets(3 * (lastrank - firstrank + 1)); // we have to keep data order
				size_t offset = 0;
				for (int r = firstrank; r <= lastrank; ++r) {
					++offset; // totalsize
					++offset; // target (always me)
					int source = rBuffer[offset++] - firstrank;
					for (int d = 0; d < 3; ++d) {
						alloffsets[3 * source + d] = rBuffer[offset++];
					}
					for (int d = 0; d < 3; ++d) {
						offset += alloffsets[3 * source + d];
					}
				}
				for (int d = 0; d < 3; ++d) {
					int sum = d * (cdistribution[info::mpi::rank - firstrank + 1] - cdistribution[info::mpi::rank - firstrank]);
					for (int r = 0; r <= lastrank - firstrank; ++r) {
						int tmp = alloffsets[3 * r + d];
						alloffsets[3 * r + d] = sum;
						sum += tmp;
					}
				}
				offset = 0;
				for (int r = firstrank; r <= lastrank; ++r) {
					++offset; // totalsize
					++offset; // target (always me)
					int source = rBuffer[offset++] - firstrank;
					int xyzsize[3] = { rBuffer[offset++], rBuffer[offset++], rBuffer[offset++] };
					for (int d = 0; d < 3; ++d) {
						if (xyzsize[d]) {
							memcpy(coordinates.data() + alloffsets[3 * source + d], rBuffer.data() + offset, sizeof(float) * xyzsize[d]);
							offset += xyzsize[d] * (sizeof(float) / sizeof(int));
						}
					}
				}
			}

			int ncoordinates = coordinates.size() / 3;
			mesh.coordinates.reserve(mesh.coordinates.size() + ncoordinates);
			for (int c = 0; c < ncoordinates; ++c) {
				mesh.coordinates.push_back(Point(coordinates[ncoordinates * 0 + c], coordinates[ncoordinates * 1 + c], coordinates[ncoordinates * 2 + c]));
			}
			mesh.nIDs.resize(mesh.coordinates.size());
			std::iota(mesh.nIDs.end() - ncoordinates, mesh.nIDs.end(), cidoffset + cdistribution[info::mpi::rank - firstrank]);
		}

		int ntotalelements = 0;
		for (; i < _elements.size() && (p + 1 == parts || _elements[i].offset < _coordinates[p + 1].offset); ++i) {
			int esize; Element::CODE code;
			switch (_elements[i].type) {
			case Elements::Type::POINT    : esize =  1; code = Element::CODE::POINT1   ; break;
			case Elements::Type::BAR2     : esize =  2; code = Element::CODE::LINE2    ; break;
			case Elements::Type::BAR3     : esize =  3; code = Element::CODE::LINE3    ; break;
			case Elements::Type::TRIA3    : esize =  3; code = Element::CODE::TRIANGLE3; break;
			case Elements::Type::TRIA6    : esize =  6; code = Element::CODE::TRIANGLE6; break;
			case Elements::Type::QUAD4    : esize =  4; code = Element::CODE::SQUARE4  ; break;
			case Elements::Type::QUAD8    : esize =  8; code = Element::CODE::SQUARE8  ; break;
			case Elements::Type::TETRA4   : esize =  4; code = Element::CODE::TETRA4   ; break;
			case Elements::Type::TETRA10  : esize = 10; code = Element::CODE::TETRA10  ; break;
			case Elements::Type::PYRAMID5 : esize =  5; code = Element::CODE::PYRAMID5 ; break;
			case Elements::Type::PYRAMID13: esize = 13; code = Element::CODE::PYRAMID13; break;
			case Elements::Type::PENTA6   : esize =  6; code = Element::CODE::PRISMA6  ; break;
			case Elements::Type::PENTA15  : esize = 15; code = Element::CODE::PRISMA15 ; break;
			case Elements::Type::HEXA8    : esize =  8; code = Element::CODE::HEXA8    ; break;
			case Elements::Type::HEXA20   : esize = 20; code = Element::CODE::HEXA20   ; break;
			case Elements::Type::NSIDED   : esize =  0; code = Element::CODE::SIZE     ; break; // not supported
			case Elements::Type::NFACED   : esize =  0; code = Element::CODE::SIZE     ; break; // not supported
			default: esize = 0; code = Element::CODE::SIZE;
			}

			size_t ebegin = std::max(_elements[i].offset, _geofile.distribution[info::mpi::rank]);
			size_t eend = std::min(_elements[i].offset + esize * sizeof(int) * _elements[i].ne, _geofile.distribution[info::mpi::rank + 1]);
			if (ebegin < eend) {
				round(ebegin, _elements[i].offset, esize * sizeof(int));
				round(eend, _elements[i].offset, esize * sizeof(int));
				std::vector<int> elements((eend - ebegin) / sizeof(int));
				memcpy(elements.data(), _geofile.begin + ebegin - _geofile.distribution[info::mpi::rank], eend - ebegin);

				int nelements = (eend - ebegin) / (sizeof(int) * esize);
				if (code != Element::CODE::POINT1) {
					ntotalelements += nelements;
					mesh.etype.resize(mesh.etype.size() + nelements, (int)code);
					mesh.esize.resize(mesh.esize.size() + nelements, esize);
					mesh.enodes.reserve(mesh.enodes.size() + elements.size());
					for (size_t n = 0; n < elements.size(); ++n) {
						mesh.enodes.push_back(elements[n] + cidoffset - 1);
					}
				} else {
					for (size_t n = 0; n < elements.size(); ++n) {
						elements[n] += cidoffset - 1;
					}
					mesh.nregions[name].assign(elements.begin(), elements.end());
				}
			} else {
				if (code == Element::CODE::POINT1) {
					mesh.nregions[name] = {};
				}
			}
		}
		einfo.push_back(eoffsetinfo{.name=name, .offset=ntotalelements, .size=ntotalelements, .totalsize=0});
		cidoffset += _coordinates[p].nn;
	}

	mesh.eIDs.resize(mesh.etype.size());

	std::vector<esint> sum, offset;
	for (size_t i = 0; i < einfo.size(); ++i) {
		offset.push_back(einfo[i].offset);
	}
	sum.resize(offset.size());
	Communication::exscan(sum, offset);
	for (size_t i = 0; i < einfo.size(); ++i) {
		einfo[i].offset = offset[i];
		einfo[i].totalsize = sum[i];
	}
	for (size_t i = 0, prevsize = 0, eidoffset = 0; i < einfo.size(); ++i) {
		if (einfo[i].totalsize) { // non-node region
			std::iota(mesh.eIDs.begin() + prevsize, mesh.eIDs.begin() + prevsize + einfo[i].size, einfo[i].offset + eidoffset);
			mesh.eregions[einfo[i].name] = std::vector<esint>(mesh.eIDs.begin() + prevsize, mesh.eIDs.begin() + prevsize + einfo[i].size);
			prevsize += einfo[i].size;
			eidoffset += einfo[i].totalsize;
		}
	}
}

void EnsightGeometry::scanASCII()
{




}
void EnsightGeometry::parseASCII(MeshBuilder &mesh)
{
	auto toLineEnd = [] (const char* &c) {
		while (*c != '\n') { c++; }
		//while (*c != ' ') { c++; }
		c++;
	};

	auto align = [] (size_t x, int y) {
		return (x - y) % 13 == 0 ? x : x + 13 - ((x - y) % 13);
	};

	auto toNextNumber = [] (const char* &c) {
		while (*c == ' ' || *c == '\n') { c++; }
		c++;
	};

	auto alignElement = [] (size_t x, int y, int z) {
		//return (x - y) % ((z * 10) + 1) == 0 ? x : x + ((z * 10) + 1) - ((x - y) % ((z * 10) + 1));
		return (x - y) % ((z * 10) + 1) == 0 ? x : x + ((z * 10) + 1) - ((x - y) % ((z * 10) + 1));
	};	
	

	auto copyStringToVec = [] (const char* &c, std::vector<char> &vector){
		size_t size = vector.size();
		vector.resize(size + 80);
		//while (*c != '\n'){
		while (*c != '\n'){	
			vector[size++] = *c++;
		}
	};
	auto removeEmptyCharacters = [] (std::string partName){
		for (size_t i = 0; i < partName.length(); i++){
			if (partName[i] == '\000' || partName[i] == ' '){
				partName.erase(partName.begin() + i, partName.end());
			}
		}
		return partName;
	};

	std::vector <esint> vecOfPointNodes;
	std::vector <esint> vecOfElementNodes;
	std::string partName;
	int esize, ecount, partID;
	Element::CODE code;
	std::vector<char> vecOfPartName;
	Communication::serialize([&] () {
		for (const char *c = _geofile.begin; c < _geofile.end; ++c) {
			if (memcmp("coordinates", c, 11) == 0) {
				toLineEnd(c);
				int n = std::atoi(c);
				toLineEnd(c);
				_coordinates.push_back(Coordinates(_geofile.distribution[info::mpi::rank] + (c - _geofile.begin), n));
			}

				///////////////////////////////
				////////////////// SCAN OBJECTS
				///////////////////////////////
			Elements::Type elementType;
			
			if (memcmp("part", c, 4) == 0){elementType = Elements::Type::POINT;
					toLineEnd(c);
					//partID = std::atoi(c);
					toLineEnd(c);
					copyStringToVec(c, _parts);

			} else if (memcmp("point", c, 5) == 0){elementType = Elements::Type::POINT;
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("quad4", c, 5) == 0){elementType = Elements::Type::QUAD4;
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("bar2", c, 4) == 0){elementType = Elements::Type::BAR2; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("tria3", c, 5) == 0){elementType = Elements::Type::TRIA3; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("tetra4", c, 6) == 0){elementType = Elements::Type::TETRA4; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("pyramid5", c, 8) == 0){elementType = Elements::Type::PYRAMID5; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("penta6", c, 6) == 0){elementType = Elements::Type::PENTA6; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("hexa8", c, 5) == 0){elementType = Elements::Type::HEXA8; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("bar3", c, 4) == 0){elementType = Elements::Type::BAR3; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("tria6", c, 5) == 0){elementType = Elements::Type::TRIA6; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("quad8", c, 5) == 0){elementType = Elements::Type::QUAD8; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("tetra10", c, 5) == 0){elementType = Elements::Type::TETRA10; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("pyramid13", c, 9) == 0){elementType = Elements::Type::PYRAMID13; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("penta15", c, 7) == 0){elementType = Elements::Type::PENTA15; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			} else if (memcmp("hexa20", c, 6) == 0){elementType = Elements::Type::HEXA20; 
				toLineEnd(c);
				ecount = std::atoi(c);
				toLineEnd(c);
				_elements.push_back(Elements(elementType, _geofile.distribution[info::mpi::rank] + (c - _geofile.begin), ecount));

			}

		}

	});

	DistributedScanner scanner;
	scanner.synchronize(_parts, _coordinates, _elements);
	
	Communication::serialize([&] () {
		printf("parts: %lu\n", _parts.size() / 80);

		printf("file: %ld - %ld\n", _geofile.distribution[info::mpi::rank], _geofile.distribution[info::mpi::rank + 1]);
		for (size_t i = 0; i < _coordinates.size(); ++i){
			printf("coo[%d]: offset=%lu, nn=%d\n", info::mpi::rank, _coordinates[i].offset, _coordinates[i].nn);
		}

		for (size_t i = 0; i < _elements.size(); ++i){
			printf("elements[%d]: offset=%lu, ne=%d, type=%d\n", info::mpi::rank, _elements[i].offset, _elements[i].ne, (int)_elements[i].type);
		}
	});

	int rank = info::mpi::rank, size = info::mpi::size, etype;
	///////////////////////////////
	////////////////////// ELEMENTS
	///////////////////////////////
	printf("start offset %d, end offset %d, MPI PROCES %d\n", _geofile.distribution[info::mpi::rank], _geofile.distribution[info::mpi::rank+1], info::mpi::rank);
	size_t parts = _parts.size() / 80;
	size_t coordinate_id_offset = 0;
	size_t element_id_offset = 0;
	for (size_t p = 0, e = 0; p < parts; ++p) {

		/////
		std::string textString(_parts.cbegin() + 80 * p, _parts.cbegin() + 80 * p + 80);
		//printf("stringFor %s\n",textString.c_str());
		textString = removeEmptyCharacters(textString);
		/////
		//printf("parse part: %s\n", textString.c_str());

		//printf("TODO: parse coordinates with index: %d, id_offset: %d\n", p, coordinate_id_offset);
		bool isNregion = false;
		std::vector<int> eids;
		for (; e < _elements.size() && (p + 1 == parts || _elements[e].offset < _coordinates[p + 1].offset); ++e){
			isNregion = false;

			if 		(_elements[e].type == Elements::Type::POINT)	{esize = 1;		etype = 0; 		isNregion = true;}
			else if (_elements[e].type == Elements::Type::BAR2)		{esize = 2;		etype = 1;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::TRIA3) 	{esize = 3;		etype = 2;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::QUAD4) 	{esize = 4;		etype = 3;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::TETRA4) 	{esize = 4;		etype = 4;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::PYRAMID5) {esize = 5;		etype = 5;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::PENTA6) 	{esize = 6;		etype = 6;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::HEXA8) 	{esize = 8;		etype = 7;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::BAR3) 	{esize = 3;		etype = 8;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::TRIA6) 	{esize = 6;		etype = 9;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::QUAD8) 	{esize = 8;		etype = 10;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::TETRA10) 	{esize = 10;	etype = 11;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::PYRAMID13){esize = 13;	etype = 12;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::PENTA15) 	{esize = 15;	etype = 13;		isNregion = false;}
			else if (_elements[e].type == Elements::Type::HEXA20) 	{esize = 20;	etype = 14;		isNregion = false;}

			size_t E = std::max(_elements[e].offset, _geofile.distribution[info::mpi::rank]);
			size_t F = std::min(_elements[e].offset + _elements[e].ne * ((esize * 10) + 1), _geofile.distribution[info::mpi::rank + 1]);
			if (E < F) {
				E = alignElement(E, _elements[e].offset, esize);
				F = alignElement(F, _elements[e].offset, esize);
				int eid_offset = (E - _elements[e].offset) / ((10 * esize) + 1);
				E -= _geofile.distribution[info::mpi::rank];
				F -= _geofile.distribution[info::mpi::rank];

				int cc = 0;
				//printf("\npush ");
				for (const char *d = _geofile.begin + E; d < _geofile.begin + F; d += 10) {
					if (isNregion == true){
						eids.push_back(atof(d) - 1 + coordinate_id_offset);
						++d;

					} else {
						mesh.enodes.push_back(atof(d) - 1 + coordinate_id_offset);
						++cc;
						//printf("%d ", mesh.enodes.back());
						if (cc % esize == 0) {
							eids.push_back(element_id_offset + eid_offset);
							mesh.eIDs.push_back(element_id_offset + eid_offset);
							++eid_offset;
							mesh.esize.push_back(esize);
							mesh.etype.push_back(etype);
							++d;
						}
						
					}
				}
			}
			element_id_offset += _elements[e].ne;

		}
		if (isNregion == true){
			mesh.nregions[textString] = eids;
		} else {
			mesh.eregions[textString] = eids;
		}
		
		// go to another part
		coordinate_id_offset += _coordinates[p].nn;
		//printf(" - - - -- - - - - - \n");
	}



	///////////////////////////////
	//////////////////////// COORDS
	/////////////////////////////// 
	std::vector<double> vectorOfXcoords;
	std::vector<double> vectorOfYcoords;
	std::vector<double> vectorOfZcoords;
	std::vector<double> vector;
	std::vector<double> vector2;
	std::vector<size_t> vecOfCoordsOffsetProcesses;
	std::vector<size_t> vecOfEndOfCoordsProcesses;
	std::vector<int> vecOfNodes;

	for (size_t i = 0; i < _coordinates.size(); ++i) {

		// najdi zacatek a konec parsovani
		size_t A = std::max(_coordinates[i].offset, _geofile.distribution[info::mpi::rank]);
		size_t B = std::min(_coordinates[i].offset + 3 * 13 * _coordinates[i].nn, _geofile.distribution[info::mpi::rank + 1]);

		if (A < B) {

			// posun na zacatek
			A = align(A, _coordinates[i].offset);
			A -= _geofile.distribution[info::mpi::rank];
			B -= _geofile.distribution[info::mpi::rank];

			// hledani cisla procesu offset a endcoords
			std::vector<size_t>::iterator low1, low2;
			low1 = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[i].offset);
			low2 = std::lower_bound(_geofile.distribution.begin(), _geofile.distribution.end(), _coordinates[i].offset + 3 * 13 * _coordinates[i].nn);
			size_t cooOffsetProcess = low1 - _geofile.distribution.begin() - 1;
			size_t cooEndProcess = low2 - _geofile.distribution.begin() - 1;
			vecOfCoordsOffsetProcesses.push_back(cooOffsetProcess);
			vecOfEndOfCoordsProcesses.push_back(cooEndProcess);



			// parsovani
			for (const char *c = _geofile.begin + A; c < _geofile.begin + B; ++c) {
				//printf("%c%c%c\n", *c, *(c  +1), *(c + 2));
				vector.push_back(atof(c));
				toLineEnd(c);
				c--;
			}
			
			int TAG = 0;

			// posilani dat
			if (rank > cooOffsetProcess && rank <= cooEndProcess) {
				MPI_Send(vector.data(), vector.size(), MPI_DOUBLE, cooOffsetProcess, TAG, MPI_COMM_WORLD);
				vector.clear();
			// prijimam
			} else if (rank == cooOffsetProcess){
				for (int r = rank + 1; r <= cooEndProcess; ++r) {
					size_t _A = std::max(_coordinates[i].offset, _geofile.distribution[r]);//r
					size_t _B = std::min(_coordinates[i].offset + 3 * 13 * _coordinates[i].nn, _geofile.distribution[r + 1]);//r+1
					_A = align(_A, _coordinates[i].offset);
					_B = align(_B, _coordinates[i].offset);
					size_t _C = (_B - _A) / 13;
					vector2.resize(_C);
					MPI_Recv(vector2.data(), vector2.size(), MPI_DOUBLE, r, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					vector.insert(vector.end(), vector2.begin(), vector2.end());

				}
			}
		}


	}
	///////////////////////////////
	///////////// ROZDELENI DLE XYZ
	/////////////////////////////// 

	for (int rankNr = 0; rankNr < info::mpi::size; rankNr++){
		if (rank == rankNr){
			//urceni poctu nodu ve vektoru
			for (int i = 0; i < _coordinates.size(); i++){
				if (_coordinates[i].offset >= _geofile.distribution[info::mpi::rank] && _coordinates[i].offset < _geofile.distribution[info::mpi::rank+1]){
					vecOfNodes.push_back(_coordinates[i].nn);
				}
			}
			int prev = 0;
			//sorting xyz
			for (int x = 0; x < vecOfNodes.size(); x++){
				for (int i = prev; i < prev + vecOfNodes[x]; i++){
					vectorOfXcoords.push_back(vector[i]);
				}
				for (int i = prev + vecOfNodes[x]; i < prev + vecOfNodes[x] * 2; i++){
					vectorOfYcoords.push_back(vector[i]);
				}	
				for (int i = prev + vecOfNodes[x] * 2; i < prev + vecOfNodes[x] * 3; i++){
					vectorOfZcoords.push_back(vector[i]);
				}
				prev = prev + vecOfNodes[x] * 3;
			}
		} 
	}

	///////////////////////////////
	////////COORDS AND nIDs TO MESH
	///////////////////////////////
	int nID = 0, prevNodes = 0;

	//spocti nody v predchozich procesech
	for (size_t i = 0; i < _coordinates.size(); i++){
		if (_coordinates[i].offset < _geofile.distribution[info::mpi::rank]){
			prevNodes = _coordinates[i].nn + prevNodes;
		}
	}
	for (size_t l = 0; l < vectorOfXcoords.size(); ++l){
		mesh.coordinates.push_back(Point(vectorOfXcoords[l], vectorOfYcoords[l], vectorOfZcoords[l]));
		mesh.nIDs.push_back(prevNodes + nID);
		nID++;
	}



	vectorOfXcoords.clear();
	vectorOfYcoords.clear();
	vectorOfZcoords.clear();
	vector.clear();
	vectorOfXcoords.resize(0);
	vectorOfYcoords.resize(0);
	vectorOfZcoords.resize(0);
	vector.resize(0);
	printf("a");
}
