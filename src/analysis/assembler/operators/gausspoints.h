
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_

#include "analysis/assembler/operator.hpp"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

template<int dimension>
struct FromNodesToGaussPoints: public Operator {
	FromNodesToGaussPoints(ParameterData &N, ParameterData &nodeData, ParameterData &gpData, int interval)
	: Operator(interval, gpData.isconst[interval], gpData.update[interval]),
	  N(N, interval, 0),
	  n(nodeData, interval),
	  gp(gpData, interval)
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

template<int dimension>
struct FromNodesToGaussPointsSimd: public Operator {
	FromNodesToGaussPointsSimd(ParameterData &N, ParameterData &nodeData, ParameterData &gpData, int interval)
	: Operator(interval, gpData.isconst[interval], gpData.update[interval]),
	  N(N, interval, 0),
	  n(nodeData, interval),
	  gp(gpData, interval)
	{

	}

	InputParameterIterator N, n;
	OutputParameterIterator gp;

	template<int nodes, int gps>
	void operator()(int gpindex)
	{	
		NtoGPSimd<nodes, dimension>(N.data + gpindex * nodes * SIMD::size, n.data, gp.data + gpindex * dimension * SIMD::size);
	}

	void operator++()
	{
		++n; ++gp;
	}

	FromNodesToGaussPointsSimd& operator+=(const int rhs)
	{
		n += rhs; gp += rhs;
		return *this;
	}
};

template <class Assembler, int dimension>
struct ElementsGaussPointsBuilder: public ElementOperatorBuilder {
	Assembler &assembler;
	ParameterData &N, &nodeData, &gpData;

	ElementsGaussPointsBuilder(Assembler &assembler, ParameterData &N, ParameterData &nodeData, ParameterData &gpData, const char *name)
	: ElementOperatorBuilder(name), assembler(assembler), N(N), nodeData(nodeData), gpData(gpData)
	{

	}

	bool build() override
	{
		gpData.addInput(nodeData);
		gpData.resize();
		assembler.addParameter(gpData);
		return true;
	}

	void apply(int interval)
	{
		iterate_elements_gps<Assembler::NGP>(FromNodesToGaussPoints<dimension>(N, nodeData, gpData, interval));
	}
};

template <class Assembler, int dimension>
struct BoundaryGaussPointsBuilder: public BoundaryOperatorBuilder {
	Assembler &assembler;
	BoundaryParameterPack &N, &nodeData, &gpData;

	BoundaryGaussPointsBuilder(Assembler &assembler, BoundaryParameterPack &N, BoundaryParameterPack &nodeData, BoundaryParameterPack &gpData, const char *name)
	: BoundaryOperatorBuilder(name), assembler(assembler), N(N), nodeData(nodeData), gpData(gpData)
	{

	}

	bool build() override
	{
		for (size_t r = 0; r < gpData.regions.size(); ++r) {
			gpData.regions[r].isset = nodeData.regions[r].isset;
			gpData.regions[r].addInput(nodeData.regions[r]);
			gpData.regions[r].resize();
			assembler.addParameter(gpData.regions[r]);
		}
		return true;
	}

	void apply(int region, int interval)
	{
		if (nodeData.regions[region].isset) {
			iterate_boundary_gps<Assembler::NGP>(FromNodesToGaussPoints<dimension>(N.regions[region], nodeData.regions[region], gpData.regions[region], interval), region);
		}
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_GAUSSPOINTS_H_ */
