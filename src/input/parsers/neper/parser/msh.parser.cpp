
#include "msh.parser.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "mesh/element.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"
#include "input/parsers/asciiparser.h"

#include <numeric>

using namespace espreso;

NeperMshMesh::Keyword::Keyword(InputFilePack &pack, const char *c)
: begin(0), end(0)
{
	rbegin = rend = info::mpi::rank;
	offset = pack.distribution[info::mpi::rank] + (c - pack.begin);
}

NeperMshMesh::KeywordWithSize::KeywordWithSize(InputFilePack &pack, const char *c, const std::string &open, const std::string &close)
: Keyword(pack, c)
{
	if (StringCompare::caseInsensitiveEq(std::string(c, c + open.size()), open)) {
		const char *_c = c;
		while (*c++ != '\n'); // $Nodes
		char *next;
		size = strtol(c, &next, 0);
		c = next;
		while (*c++ != '\n'); // go to the first data line
		begin = offset + (c - _c);
	} else if (StringCompare::caseInsensitiveEq(std::string(c, c + close.size()), close)) {
		const char *_c = c;
		while (*c != '\n') { --c; }
		end = offset + (c - _c);
	} else {
		eslog::error("MESIO internal error: cannot parse 'NEPER.msh format.\n");
	}
}

NeperMshMesh::Format::Format(InputFilePack &pack, const char *c)
: Keyword(pack, c)
{
	if (StringCompare::caseInsensitiveEq(std::string(c, c + 11), "$MeshFormat")) {
		const char *_c = c;
		while (*c++ != '\n'); // $MeshFormat
		begin = offset + (c - _c);
		char *next;
		version = strtod(c, &next);
		c = next;
		file_type = strtol(c, &next, 0);
		c = next;
		data_size = strtol(c, &next, 0);
	} else if (StringCompare::caseInsensitiveEq(std::string(c, c + 14), "$EndMeshFormat")) {
		const char *_c = c;
		while (*c != '\n') { --c; }
		end = offset + (c - _c);
	} else {
		eslog::error("MESIO internal error: cannot parse 'NEPER.msh format.\n");
	}
}

NeperMshMesh::NeperMshMesh(InputFilePack &meshfile)
: _meshfile(meshfile)
{

}

void NeperMshMesh::parse(MeshBuilder &mesh)
{
	profiler::syncstart("neper_parse");
	profiler::syncparam("size", _meshfile.end - _meshfile.begin);

	scan();
	parseASCII(mesh);

	profiler::syncend("neper_parse");
}

template<typename T>
void add(DistributedScanner &scanner, InputFilePack &meshfile, std::vector<T> &data)
{
	scanner.add({ T::open, T::close }, [&] (const char *c) { data.push_back(T(meshfile, c)); });
}

template<typename T>
void merge(std::vector<T> &keywords)
{
	if (keywords.size() == 2) {
		keywords.front().end = keywords.back().end;
		keywords.front().rend = keywords.back().rend;
		keywords.pop_back();
	} else if (!keywords.empty()) {
		eslog::error("MESIO internal error: cannot parse 'NEPER.msh format.\n");
	}
}

void NeperMshMesh::scan()
{
	DistributedScanner scanner;

	add<Format>(scanner, _meshfile, _format);
	add<Nodes>(scanner, _meshfile, _nodes);
	add<Elements>(scanner, _meshfile, _elements);
	add<Periodicity>(scanner, _meshfile, _periodicity);
	add<NSets>(scanner, _meshfile, _nsets);
	add<Fasets>(scanner, _meshfile, _fasets);
	add<NodePartitions>(scanner, _meshfile, _nodePartitions);
	add<PhysicalNames>(scanner, _meshfile, _physicalNames);
	add<ElsetOrientations>(scanner, _meshfile, _elsetOrientations);
	add<ElementOrientations>(scanner, _meshfile, _elementOrientations);
	add<Groups>(scanner, _meshfile, _groups);

	scanner.align(_meshfile, " \n");
	scanner.scanlines(_meshfile);

	scanner.synchronize(_format, _nodes, _elements, _periodicity, _nsets, _fasets, _nodePartitions, _physicalNames, _elsetOrientations, _elementOrientations, _groups);

	merge(_format);
	merge(_nodes);
	merge(_elements);
	merge(_periodicity);
	merge(_nsets);
	merge(_fasets);
	merge(_nodePartitions);
	merge(_physicalNames);
	merge(_elsetOrientations);
	merge(_elementOrientations);
	merge(_groups);
}

struct Boundary {
	char name[NEPER_MAX_NAME_SIZE] = { 0 };
	int rank = info::mpi::rank;
};

void NeperMshMesh::parseASCII(MeshBuilder &mesh)
{
	// nodes
	ASCIIParser::parse(mesh.nIDs, mesh.coordinates, _meshfile, _nodes.front().begin, _nodes.front().end);

	// elements
	std::vector<esint> edata;
	std::vector<std::vector<esint> > eregs;
	ASCIIParser::parse(edata, _meshfile, _elements.front().begin, _elements.front().end);
	if (_elements.front().rend != _elements.front().rbegin) {
		mesh.eIDs.reserve(_elements.front().size / (_elements.front().rend - _elements.front().rbegin));
	} else {
		mesh.eIDs.reserve(_elements.front().size);
	}
	mesh.etype.reserve(mesh.eIDs.capacity());
	mesh.esize.reserve(mesh.eIDs.capacity());
	mesh.enodes.reserve(mesh.eIDs.capacity() * 8);
	for (size_t i = 0; i < edata.size(); ) {
		esint id, type, dim, tags, enodes = 0;
		size_t tag0, tag1;
		id = edata[i++];
		switch (edata[i++]) {
		case 15: dim = 0; type = (int)Element::CODE::POINT1;    enodes =  1; break;
		case 1:  dim = 1; type = (int)Element::CODE::LINE2;     enodes =  2; break;
		case 8:  dim = 1; type = (int)Element::CODE::LINE3;     enodes =  3; break;
		case 2:  dim = 2; type = (int)Element::CODE::TRIANGLE3; enodes =  3; break;
		case 3:  dim = 2; type = (int)Element::CODE::SQUARE4;   enodes =  4; break;
		case 9:  dim = 2; type = (int)Element::CODE::TRIANGLE6; enodes =  6; break;
		case 16: dim = 2; type = (int)Element::CODE::SQUARE8;   enodes =  8; break;
		case 4:  dim = 3; type = (int)Element::CODE::TETRA4;    enodes =  4; break;
		case 5:  dim = 3; type = (int)Element::CODE::HEXA8;     enodes =  8; break;
		case 11: dim = 3; type = (int)Element::CODE::TETRA10;   enodes = 10; break;
		case 17: dim = 3; type = (int)Element::CODE::HEXA20;    enodes = 20; break;
		case 6:  dim = 3; type = (int)Element::CODE::PRISMA6;   enodes =  6; break;
		case 18: dim = 3; type = (int)Element::CODE::PRISMA15;  enodes = 15; break;
		case 10: eslog::error("MESIO internal error: unsupported element SQUARE with 9 nodes."); break;
		}
		switch (tags = edata[i++]) { // tags
		case 0: break;
		case 1: tag0 = edata[i++]; break;
		case 2: tag0 = edata[i++]; tag1 = edata[i++]; break;
		case 3: tag0 = edata[i++]; tag1 = edata[i++]; i++; break; // tag2 = partition -> skipped
		default:
			eslog::error("MESIO internal error: unsupported number of tags (more than 3).");
		}
		if (dim == 3) {
			mesh.eIDs.push_back(id);
			mesh.etype.push_back(type);
			mesh.esize.push_back(enodes);
			for (esint n = 0; n < enodes; ++n) {
				mesh.enodes.push_back(edata[i++]);
			}
			if (tag0 == tag1) {
				if (tag0 >= eregs.size()) {
					eregs.resize(tag0);
				}
				eregs[tag0 - 1].push_back(id);
			} else {
				eslog::error("MESIO internal error: tag0 != tag1.");
			}
		} else {
			i += enodes;
		}
	}

	std::vector<Boundary> boundaries;

	// sets
	std::vector<esint> nsets;
	ASCIIParser::parseWithStrings(nsets, _meshfile, _nsets.front().begin, _nsets.front().end);
	for (size_t i = 0; i < nsets.size(); ++i) {
		if (nsets[i] < 0) {
			boundaries.push_back(Boundary());

			const char *c = _meshfile.begin + (-nsets[i]);
			char *to = boundaries.back().name;
			while (*c != '\n') { *to++ = *c++; }
		}
	}

	std::vector<esint> polyid;
	std::vector<Point> orientation;
	ASCIIParser::parse(polyid, orientation, _meshfile, _elsetOrientations.front().begin, _elsetOrientations.front().end);

	esint maxpoly = eregs.size();
	Communication::allReduce(&maxpoly, NULL, 1, MPITools::getType(maxpoly).mpitype, MPI_MAX);
	Communication::allGatherUnknownSize(boundaries);
	Communication::allGatherUnknownSize(polyid);
	Communication::allGatherUnknownSize(orientation);

	size_t current = 0;
	for (size_t i = 0; i < boundaries.size(); ++i) {
		mesh.nregions[boundaries[i].name];
		if (boundaries[i].rank < info::mpi::rank) {
			current = i;
		}
	}

	size_t begin = 0;
	for (size_t i = 0; i < nsets.size(); ++i) {
		if (nsets[i] < 0) {
			if (i) {
				mesh.nregions[boundaries[current++].name].assign(nsets.begin() + begin, nsets.begin() + i);
			}
			++i; // skip number of nodes
			begin = i + 1;
		}
	}
	mesh.nregions[boundaries[current].name].assign(nsets.begin() + begin, nsets.end());

	eregs.resize(maxpoly);
	for (esint i = 0; i < maxpoly; ++i) {
		mesh.eregions["poly" + std::to_string(i + 1)].swap(eregs[i]);
	}

	for (size_t i = 0; i < polyid.size(); ++i) {
		mesh.orientation["poly" + std::to_string(polyid[i])] = orientation[i];
	}
}


