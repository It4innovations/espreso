
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
	}

	CoordinatesToElementNodes& operator+=(const int rhs)
	{
		procNodes += rhs;
		return *this;
	}
};

struct Coordinates2DToElementNodes: CoordinatesToElementNodes {
	using CoordinatesToElementNodes::CoordinatesToElementNodes;

	void operator()()
	{
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
			ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
			ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
		}
	}
};

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
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_ */
