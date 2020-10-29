
#ifndef SRC_MESH_ELEMENT_H_
#define SRC_MESH_ELEMENT_H_

#include <cstddef>
#include <vector>

namespace espreso {

struct BaseFunctions;
class MatrixDense;
template <typename TEBoundaries, typename TEData> class serializededata;
template <typename TEData> class edata;

struct Element {

	enum class TYPE: int {
		POINT  = 0,
		LINE   = 1,
		PLANE  = 2,
		VOLUME = 3,
	};

	enum class CODE: int {
		POINT1, // 0

		// without mid-points
		LINE2, // 1

		TRIANGLE3, // 2
		SQUARE4, // 3

		TETRA4, // 4
		PYRAMID5, // 5
		PRISMA6, // 6
		HEXA8, // 7

		// with mid-points
		LINE3, // 8

		TRIANGLE6, // 9
		SQUARE8, // 10

		TETRA10, // 11
		PYRAMID13, // 12
		PRISMA15, // 13
		HEXA20, // 14

		NOT_SUPPORTED,

		// number of element types
		SIZE
	};

	TYPE type;
	CODE code;
	int nodes;
	int coarseNodes;
	int nCommonFace;
	int nCommonEdge;

	std::vector<MatrixDense> *N;
	std::vector<MatrixDense> *NN;
	std::vector<MatrixDense> *NNN;
	std::vector<MatrixDense> *dN;
	std::vector<double> *weighFactor;

	std::vector<MatrixDense> *nN;
	std::vector<MatrixDense> *ndN;

	serializededata<int, int> *faces;
	serializededata<int, int> *edges;

	serializededata<int, Element*> *facepointers;
	serializededata<int, Element*> *edgepointers;

	serializededata<int, int> *triangles;
	std::vector<int> *polygon;

	Element()
	: type(TYPE::POINT), code(CODE::POINT1), nodes(1), coarseNodes(1), nCommonFace(1), nCommonEdge(1),
	  N(NULL), NN(NULL), NNN(NULL), dN(NULL), weighFactor(NULL), nN(NULL), ndN(NULL),
	  faces(NULL), edges(NULL), facepointers(NULL), edgepointers(NULL),
	  triangles(NULL), polygon(NULL) {}

	template<CODE> void init();

	virtual ~Element();

	int getIndex(edata<esint> &enodes, serializededata<int, int> *subindices, serializededata<int, Element*> *subpointers, edata<esint> &subnodes);
};

}

#endif /* SRC_MESH_ELEMENT_H_ */
