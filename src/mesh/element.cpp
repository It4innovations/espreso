
#include "element.h"

#include "basis/containers/serializededata.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

using namespace espreso;

Element::~Element()
{
	if (faces != NULL) { delete faces; }
	if (edges != NULL) { delete edges; }
	if (facepointers != NULL) { delete facepointers; }
	if (edgepointers != NULL) { delete edgepointers; }
	if (triangles != NULL) { delete triangles; }
}

int Element::getIndex(edata<esint> &enodes, serializededata<int, int> *subindices, serializededata<int, Element*> *subpointers, edata<esint> &subnodes)
{
	int index = 0;
	for (auto subi = subindices->cbegin(); subi != subindices->cend(); ++subi, ++index) {
		int nsize = subpointers->datatarray()[index]->coarseNodes;
		for (int n = 0; n < nsize; n++) {
			if (subnodes[n] == enodes[*subi->begin()]) { // find the same node
				if (subnodes[(n + 1) % nsize] == enodes[*(subi->begin() + 1)]) { // check the direction
					int nn;
					for (nn = 2; nn < nsize; nn++) {
						if (subnodes[(n + nn) % nsize] != enodes[*(subi->begin() + nn)]) {
							break;
						}
					}
					if (nn == nsize) {
						return index;
					}
				}
				if (subnodes[(n + nsize - 1) % nsize] == enodes[*(subi->begin() + 1)]) { // check the direction
					int nn;
					for (nn = 2; nn < nsize; nn++) {
						if (subnodes[(n + nsize - nn) % nsize] != enodes[*(subi->begin() + nn)]) {
							break;
						}
					}
					if (nn == nsize) {
						return index;
					}
				}
				break;
			}
		}
	}

	return -1;
}

namespace espreso {

template<> void Element::init<Element::CODE::POINT1>()
{
	type = Element::TYPE::POINT;
	code = Element::CODE::POINT1;
	nodes = 1;
	coarseNodes = 1;
	nCommonFace = 1;
	nCommonEdge = 1;
}

template<> void Element::init<Element::CODE::LINE2>()
{
	type = Element::TYPE::LINE;
	code = Element::CODE::LINE2;
	nodes = 2;
	coarseNodes = 2;
	nCommonFace = 1;
	nCommonEdge = 1;
}

template<> void Element::init<Element::CODE::LINE3>()
{
	type = Element::TYPE::LINE;
	code = Element::CODE::LINE3;
	nodes = 3;
	coarseNodes = 2;
	nCommonFace = 1;
	nCommonEdge = 1;
}

template<> void Element::init<Element::CODE::TRIANGLE3>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::TRIANGLE3;
	nodes = 3;
	coarseNodes = 3;
	nCommonFace = 2;
	nCommonEdge = 1;

	std::vector<Element*> epointers(3, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);

	std::vector<int> data = {
		0, 1,
		1, 2,
		2, 0
	};

	std::vector<int> tringles = {
		0, 1, 2
	};

	edges = new serializededata<int, int>(2, data);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(2, data);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);
}

template<> void Element::init<Element::CODE::TRIANGLE6>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::TRIANGLE6;
	nodes = 6;
	coarseNodes = 3;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> epointers(3, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);

	std::vector<int> data = {
		0, 1, 3,
		1, 2, 4,
		2, 0, 5
	};

	std::vector<int> tringles = {
		0, 3, 5,
		3, 1, 4,
		4, 2, 5,
		5, 0, 3,
		3, 4, 5
	};

	edges = new serializededata<int, int>(3, data);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(3, data);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);
}

template<> void Element::init<Element::CODE::SQUARE4>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::SQUARE4;
	nodes = 4;
	coarseNodes = 4;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> epointers(4, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);

	std::vector<int> lines = {
		0, 1,
		1, 2,
		2, 3,
		3, 0
	};

	std::vector<int> tringles = {
		0, 1, 2,
		0, 2, 3
	};

	edges = new serializededata<int, int>(2, lines);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(2, lines);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);
}

template<> void Element::init<Element::CODE::SQUARE8>()
{
	type = Element::TYPE::PLANE;
	code = Element::CODE::SQUARE8;
	nodes = 8;
	coarseNodes = 4;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> epointers(4, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);

	std::vector<int> data = {
		0, 1, 4,
		1, 2, 5,
		2, 3, 6,
		3, 0, 7
	};

	std::vector<int> tringles = {
		0, 4, 7,
		4, 1, 5,
		5, 2, 6,
		6, 3, 7,
		4, 5, 6,
		4, 6, 7
	};

	edges = new serializededata<int, int>(3, data);
	edgepointers = new serializededata<int, Element*>(1, epointers);
	faces = new serializededata<int, int>(3, data);
	facepointers = new serializededata<int, Element*>(1, epointers);
	triangles = new serializededata<int, int>(3, tringles);
}

template<> void Element::init<Element::CODE::TETRA4>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::TETRA4;
	nodes = 4;
	coarseNodes = 4;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers(4, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE3)]);

	std::vector<int> fpoints = {
		0, 1, 3,
		1, 2, 3,
		2, 0, 3,
		2, 1, 0
	};

	faces = new serializededata<int, int>(3, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(6, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);
	std::vector<int> epoints = {
		0, 1,
		1, 2,
		2, 0,
		0, 3,
		1, 3,
		2, 3
	};

	edges = new serializededata<int, int>(2, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);
}

template<> void Element::init<Element::CODE::TETRA10>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::TETRA10;
	nodes = 10;
	coarseNodes = 4;
	nCommonFace = 4;
	nCommonEdge = 3;

	std::vector<Element*> fpointers(4, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE6)]);

	std::vector<int> fpoints = {
		0, 1, 3, 4, 8, 7,
		1, 2, 3, 5, 9, 8,
		2, 0, 3, 6, 7, 9,
		2, 1, 0, 5, 4, 6
	};

	faces = new serializededata<int, int>(6, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(6, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);

	std::vector<int> epoints = {
		0, 1, 4,
		1, 2, 5,
		2, 0, 6,
		0, 3, 7,
		1, 3, 8,
		2, 3, 9,
	};

	edges = new serializededata<int, int>(3, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);
}

template<> void Element::init<Element::CODE::PYRAMID5>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PYRAMID5;
	nodes = 5;
	coarseNodes = 5;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers;
	fpointers.resize(1, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE4)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE3)]);

	std::vector<int> fdist = { 0, 4, 7, 10, 13, 16 };
	std::vector<int> fpoints = {
		3, 2, 1, 0,
		0, 1, 4,
		1, 2, 4,
		2, 3, 4,
		3, 0, 4
	};

	faces = new serializededata<int, int>(fdist, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(8, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);
	std::vector<int> epoints = {
		0, 1,
		1, 2,
		2, 3,
		3, 0,
		0, 4,
		1, 4,
		2, 4,
		3, 4,
	};

	edges = new serializededata<int, int>(2, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);
}

template<> void Element::init<Element::CODE::PYRAMID13>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PYRAMID13;
	nodes = 13;
	coarseNodes = 5;
	nCommonFace = 4;
	nCommonEdge = 3;

	std::vector<Element*> fpointers;
	fpointers.resize(1, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE8)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE6)]);

	std::vector<int> fdist = { 0, 8, 14, 20, 26, 32 };
	std::vector<int> fpoints = {
		3, 2, 1, 0,  7,  6, 5, 8,
		0, 1, 4, 5, 10,  9,
		1, 2, 4, 6, 11, 10,
		2, 3, 4, 7, 12, 11,
		3, 0, 4, 8,  9, 12
	};

	faces = new serializededata<int, int>(fdist, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(8, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);

	std::vector<int> epoints = {
		0, 1,  5,
		1, 2,  6,
		2, 3,  7,
		3, 0,  8,
		0, 4,  9,
		1, 4, 10,
		2, 4, 11,
		3, 4, 12,
	};

	edges = new serializededata<int, int>(3, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);
}

template<> void Element::init<Element::CODE::PRISMA6>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PRISMA6;
	nodes = 6;
	coarseNodes = 6;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers;
	fpointers.resize(3, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE4)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE3)]);

	std::vector<int> fdist = { 0, 4, 8, 12, 15, 18 };
	std::vector<int> fpoints = {
		0, 1, 4, 3,
		1, 2, 5, 4,
		2, 0, 3, 5,
		2, 1, 0,
		3, 4, 5,
	};

	faces = new serializededata<int, int>(fdist, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(9, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);
	std::vector<int> epoints = {
		0, 1,
		1, 2,
		2, 0,
		3, 4,
		4, 5,
		5, 3,
		0, 3,
		1, 4,
		2, 5,
	};

	edges = new serializededata<int, int>(2, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);
}

template<> void Element::init<Element::CODE::PRISMA15>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::PRISMA15;
	nodes = 15;
	coarseNodes = 6;
	nCommonFace = 4;
	nCommonEdge = 3;

	std::vector<Element*> fpointers;
	fpointers.resize(3, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE8)]);
	fpointers.resize(5, &Mesh::edata[static_cast<int>(Element::CODE::TRIANGLE6)]);

	std::vector<int> fdist = { 0, 8, 16, 24, 30, 36 };
	std::vector<int> fpoints = {
		0, 1, 4, 3,  6, 13,  9, 12,
		1, 2, 5, 4,  7, 14, 10, 13,
		2, 0, 3, 5,  8, 12, 11, 14,
		2, 1, 0, 7,  6,  8,
		3, 4, 5, 9, 10, 11
	};

	faces = new serializededata<int, int>(fdist, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(9, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);
	std::vector<int> epoints = {
		0, 1,  6,
		1, 2,  7,
		2, 0,  8,
		3, 4,  9,
		4, 5, 10,
		5, 3, 11,
		0, 3, 12,
		1, 4, 13,
		2, 5, 14,
	};

	edges = new serializededata<int, int>(3, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);
}

template<> void Element::init<Element::CODE::HEXA8>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::HEXA8;
	nodes = 8;
	coarseNodes = 8;
	nCommonFace = 3;
	nCommonEdge = 2;

	std::vector<Element*> fpointers(6, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE4)]);
	std::vector<int> fpoints = {
		0, 1, 5, 4,
		3, 2, 1, 0,
		4, 5, 6, 7,
		7, 6, 2, 3,
		1, 2, 6, 5,
		3, 0, 4, 7
	};

	faces = new serializededata<int, int>(4, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(12, &Mesh::edata[static_cast<int>(Element::CODE::LINE2)]);
	std::vector<int> epoints = {
		0, 1,
		1, 2,
		2, 3,
		3, 0,
		4, 5,
		5, 6,
		6, 7,
		7, 4,
		0, 4,
		1, 5,
		2, 6,
		3, 7,
	};

	edges = new serializededata<int, int>(2, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);

}

template<> void Element::init<Element::CODE::HEXA20>()
{
	type = Element::TYPE::VOLUME;
	code = Element::CODE::HEXA20;
	nodes = 20;
	coarseNodes = 8;
	nCommonFace = 4;
	nCommonEdge = 3;

	std::vector<Element*> fpointers(6, &Mesh::edata[static_cast<int>(Element::CODE::SQUARE8)]);

	std::vector<int> fpoints = {
		0, 1, 5, 4,  8, 17, 12, 16,
		3, 2, 1, 0, 10,  9,  8, 11,
		4, 5, 6, 7, 12, 13, 14, 15,
		7, 6, 2, 3, 14, 18, 10, 19,
		1, 2, 6, 5,  9, 18, 13, 17,
		3, 0, 4, 7, 11, 16, 15, 19
	};

	faces = new serializededata<int, int>(8, fpoints);
	facepointers = new serializededata<int, Element*>(1, fpointers);

	std::vector<Element*> epointers(12, &Mesh::edata[static_cast<int>(Element::CODE::LINE3)]);

	std::vector<int> epoints = {
		0, 1,  8,
		1, 2,  9,
		2, 3, 10,
		3, 0, 11,
		4, 5, 12,
		5, 6, 13,
		6, 7, 14,
		7, 4, 15,
		0, 4, 16,
		1, 5, 17,
		2, 6, 18,
		3, 7, 19,
	};

	edges = new serializededata<int, int>(3, epoints);
	edgepointers = new serializededata<int, Element*>(1, epointers);
}

}
