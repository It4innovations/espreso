
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"

namespace espreso {

struct CoordinatesToElementNodes: public ActionOperator {
	CoordinatesToElementNodes(int interval, serializededata<esint, esint>::const_iterator procNodes, ParameterData &ncoordinates)
	: procNodes(procNodes),
	  ncoordinates(ncoordinates, interval, ndim)
	{

	}

	serializededata<esint, esint>::const_iterator procNodes;
	OutputParameterIterator ncoordinates;

	void operator++()
	{
		++procNodes;
	}

	void move(int n)
	{
		procNodes += n;
		ncoordinates += n * procNodes->size();
	}

	CoordinatesToElementNodes& operator+=(const int rhs)
	{
		procNodes += rhs;
		return *this;
	}
};

template <size_t nodes, size_t gps>
struct Coordinates2DToElementNodes: CoordinatesToElementNodes {
	using CoordinatesToElementNodes::CoordinatesToElementNodes;

	//   element 0    element 1
	// [xy xy xy xy][xy xy xy xy]..
	void operator()()
	{
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
			ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
			ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
		}
	}

	//   element 0, 1           element 2, 3
	// [xxyy xxyy xxyy xxyy][xxyy xxyy xxyy xxyy]
	void simd()
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				ncoordinates[(2 * n + 0) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)].x;
				ncoordinates[(2 * n + 1) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)].y;
			}
		}
		ncoordinates += SIMD::size * nodes;
	}

	void peel(size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				ncoordinates[(2 * n + 0) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)].x;
				ncoordinates[(2 * n + 1) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)].y;
			}
		}
		ncoordinates += size * nodes;
	}
};

template <size_t nodes, size_t gps>
struct Coordinates3DToElementNodes: CoordinatesToElementNodes {
	using CoordinatesToElementNodes::CoordinatesToElementNodes;

	void operator()()
	{
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
			ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
			ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
			ncoordinates[2] = info::mesh->nodes->coordinates->datatarray()[*n].z;
		}
	}

	void simd()
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				ncoordinates[(2 * n + 0) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)].x;
				ncoordinates[(3 * n + 1) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)].y;
				ncoordinates[(3 * n + 2) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)].z;
			}
		}
		ncoordinates += SIMD::size * nodes;
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_ */
