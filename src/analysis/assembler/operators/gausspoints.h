
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

template<size_t nodes, size_t gps, size_t dimension>
struct FromNodesToGaussPoints: public ActionOperator {
	FromNodesToGaussPoints(int interval, ParameterData &N, ParameterData &nodeData, ParameterData &gpData)
	: ActionOperator(interval, gpData.isconst[interval], gpData.update[interval]),
	  N(N, interval, 0),
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

	void reset()
	{

	}
};

//template <int dimension>
//struct ElementsGaussPointsBuilder: public ElementOperatorBuilder {
//	ParameterData &N, &nodeData, &gpData;
//
//	ElementsGaussPointsBuilder(AX_HeatTransfer &module, ParameterData &N, ParameterData &nodeData, ParameterData &gpData, const char *name)
//	: ElementOperatorBuilder(name), N(N), nodeData(nodeData), gpData(gpData)
//	{
//		gpData.addInput(nodeData);
//		gpData.resize();
//		module.addParameter(gpData);
//	}
//
//	void apply(int interval)
//	{
////		iterate_elements_gps<Assembler::NGP>(FromNodesToGaussPoints<dimension>(N, nodeData, gpData, interval));
//	}
//};

//template <int dimension>
//struct BoundaryGaussPointsBuilder: public BoundaryOperatorBuilder {
//	BoundaryParameterPack &N, &nodeData, &gpData;
//
//	BoundaryGaussPointsBuilder(AX_HeatTransfer &module, BoundaryParameterPack &N, BoundaryParameterPack &nodeData, BoundaryParameterPack &gpData, const char *name)
//	: BoundaryOperatorBuilder(name), N(N), nodeData(nodeData), gpData(gpData)
//	{
//		for (size_t r = 0; r < gpData.regions.size(); ++r) {
//			gpData.regions[r].isset = nodeData.regions[r].isset;
//			gpData.regions[r].addInput(nodeData.regions[r]);
//			gpData.regions[r].resize();
//			module.addParameter(gpData.regions[r]);
//		}
//	}
//
//	void apply(int region, int interval)
//	{
////		if (nodeData.regions[region].isset) {
////			iterate_boundary_gps<Assembler::NGP>(FromNodesToGaussPoints<dimension>(N.regions[region], nodeData.regions[region], gpData.regions[region], interval), region);
////		}
//	}
//};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_ */
