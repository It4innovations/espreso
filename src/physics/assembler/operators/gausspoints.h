
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/math.hpp"

namespace espreso {

template<int dimension>
struct FromNodesToGaussPoints: public Operator {
	GET_NAME(FromNodesToGaussPoints)

	FromNodesToGaussPoints(ParameterData &N, ParameterData &nodeData, ParameterData &gpData, int interval)
	: Operator(interval, gpData.isconst[interval], Link(interval).inputs(N, nodeData).outputs(gpData)),
	  N(N, interval, 0),
	  n(nodeData, interval, dimension * enodes),
	  gp(gpData, interval, dimension * egps)
	{

	}

	InputParameterIterator N, n;
	OutputParameterIterator gp;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{
		NtoGP<nodes, dimension>(N.data + gpindex * nodes, n.data, gp.data + gpindex * dimension);
	}

	void operator++()
	{
		++n; ++gp;
	}
};

template <class Kernel, int dimension>
struct ElementsGaussPointsBuilder: public ElementOperatorBuilder {
	GET_NAME(ElementsGaussPointsBuilder)

	ParameterData &N, &nodeData, &gpData;

	ElementsGaussPointsBuilder(ParameterData &N, ParameterData &nodeData, ParameterData &gpData)
	: N(N), nodeData(nodeData), gpData(gpData)
	{
		gpData.addInputs(nodeData);
	}

	void apply(int interval)
	{
		iterate_elements_gps<Kernel>(FromNodesToGaussPoints<dimension>(N, nodeData, gpData, interval));
	}
};

template <class Kernel, int dimension>
struct BoundaryGaussPointsBuilder: public BoundaryOperatorBuilder {
	GET_NAME(BoundaryGaussPointsBuilder)

	BoundaryParameterPack &N, &nodeData, &gpData;

	BoundaryGaussPointsBuilder(BoundaryParameterPack &N, BoundaryParameterPack &nodeData, BoundaryParameterPack &gpData)
	: N(N), nodeData(nodeData), gpData(gpData)
	{
		for (size_t r = 0; r < gpData.regions.size(); ++r) {
			gpData.regions[r].isset = nodeData.regions[r].isset;
			gpData.regions[r].addInputs(nodeData.regions[r]);
		}
	}

	void apply(int region, int interval)
	{
		if (nodeData.regions[region].isset) {
			iterate_boundary_gps<Kernel>(FromNodesToGaussPoints<dimension>(N.regions[region], nodeData.regions[region], gpData.regions[region], interval), region);
		}
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_ */
