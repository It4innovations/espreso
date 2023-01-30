
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_

#include "analysis/assembler/operator.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics>
struct CoordinatesToElementNodes: ActionOperator, Physics {
	serializededata<esint, esint>::const_iterator procNodes;

	CoordinatesToElementNodes(int interval, serializededata<esint, esint>::const_iterator procNodes)
	: procNodes(procNodes)
	{

	}

	void move(int n)
	{
		procNodes += n;
	}

	//   element 0    element 1
	// [xy xy xy xy][xy xy xy xy]..
	void sisd(typename Physics::Element &element)
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t d = 0; d < ndim; ++d) {
				element.coords[ndim * n + d] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
			}
		}
		++procNodes;
	}

	//   element 0, 1           element 2, 3
	// [xxyy xxyy xxyy xxyy][xxyy xxyy xxyy xxyy]
	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.coords[(ndim * n + d) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
				}
			}
		}
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.coords[(ndim * n + d) * SIMD::size + s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
				}
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct CoordinatesToElementNodes2: ActionOperator, Physics {
	serializededata<esint, esint>::const_iterator procNodes;

	CoordinatesToElementNodes2(int interval, serializededata<esint, esint>::const_iterator procNodes)
	: procNodes(procNodes)
	{
		procNodes += info::mesh->elements->eintervals[interval].begin;
		isconst = false;
	}

	void move(int n)
	{
		procNodes += n;
	}

	//   element 0    element 1
	// [xy xy xy xy][xy xy xy xy]..
	void sisd(typename Physics::Element &element)
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t d = 0; d < ndim; ++d) {
				element.coords[n][d][0] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
			}
		}
		++procNodes;
	}

	//   element 0, 1           element 2, 3
	// [xxyy xxyy xxyy xxyy][xxyy xxyy xxyy xxyy]
	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0; s < SIMD::size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
				}
			}
		}
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0; s < size; ++s, ++procNodes) {
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.coords[n][d][s] = info::mesh->nodes->coordinates->datatarray()[procNodes->at(n)][d];
				}
			}
		}
	}
};


}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_ */
