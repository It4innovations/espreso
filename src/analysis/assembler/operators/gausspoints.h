
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_

#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

template<size_t nodes, size_t gps, size_t dimension>
struct FromNodesToGaussPoints: public ActionOperator {
	FromNodesToGaussPoints(int interval, const ParameterData &N, const ParameterData &nodeData, ParameterData &gpData)
	: N(N, interval, 0),
	  n(nodeData, interval),
	  gp(gpData, interval)
	{

	}

	InputParameterIterator N, n;
	OutputParameterIterator gp;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			NtoGP<nodes, dimension>(N.data + gpindex * nodes, n.data, gp.data + gpindex * dimension);
		}
	}

	void operator++()
	{
		++n; ++gp;
	}

	void move(int n)
	{
		this->n += n; gp += n;
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_ */
