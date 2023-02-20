
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_

#include "analysis/assembler/operator.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "math/simd/simd.h"
#include "mesh/store/nodestore.h"

namespace espreso {

struct CopyCoordinates: ActionOperator {
	serializededata<esint, esint>::const_iterator procNodes;

	CopyCoordinates(size_t interval, serializededata<esint, esint>::const_iterator procNodes)
	: procNodes(procNodes)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	CopyCoordinates(size_t region, size_t interval, serializededata<esint, esint>::const_iterator procNodes)
	: procNodes(procNodes)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::SOLUTION;
	}

	void move(int n)
	{
		procNodes += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct CoordinatesToElementNodes: CopyCoordinates, Physics {
	using CopyCoordinates::CopyCoordinates;

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

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct CoordinatesToElementNodesAndGPs: CopyCoordinates, Physics {
	using CopyCoordinates::CopyCoordinates;

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
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t d = 0; d < ndim; ++d) {
				element.gpcoords[gp][d] = zeros();
			}
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.gpcoords[gp][d] = element.gpcoords[gp][d] + load1(element.N[gp][n]) * element.coords[n][d];
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
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t d = 0; d < ndim; ++d) {
				element.gpcoords[gp][d] = zeros();
			}
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t d = 0; d < ndim; ++d) {
					element.gpcoords[gp][d] = element.gpcoords[gp][d] + load1(element.N[gp][n]) * element.coords[n][d];
				}
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_ */
