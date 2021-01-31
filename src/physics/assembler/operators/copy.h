
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

#include <cstdio>
#include <map>

namespace espreso {

class NodeData;
class ECFExpression;

template <class TParent>
struct CopyParameters: public TParent {
	const ParameterData &from;
	ParameterData &to;

	CopyParameters(const ParameterData &from, ParameterData &to)
	: from(from), to(to)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void apply(int interval)
	{
		if (to.version[interval] < from.version[interval]) {
			if (Operator::print) printf("\tOP::PARAMETER::%d::COPY\n", interval);
			to.version[interval] = from.version[interval];
			if (from.isconst[interval]) {
				auto indata = from.data->begin() + interval;
				auto outdata = to.data->begin() + interval;
				for (size_t i = 0; i < outdata->size() / indata->size(); ++i) {
					memcpy(outdata->data() + i * indata->size(), indata->data(), sizeof(double) * indata->size());
				}
			} else {
				memcpy((to.data->begin() + interval)->data(), (from.data->begin() + interval)->data(), sizeof(double) * (from.data->begin() + interval)->size());
			}
		} else {
			if (Operator::print > 1) printf("\tOP::PARAMETER::COPY::%d::SKIPPED\n", interval);
		}
	}
};

struct CopyElementParameters: public CopyParameters<ElementOperatorBuilder> {
	GET_NAME(CopyElementParameters)
	using CopyParameters<ElementOperatorBuilder>::CopyParameters;
};
struct CopyBoundaryParameters: public CopyParameters<BoundaryOperatorBuilder> {
	GET_NAME(CopyBoundaryParameters)
	using CopyParameters<BoundaryOperatorBuilder>::CopyParameters;
};

struct CopyNodesToElementsNodes: public OperatorBuilder {
	GET_NAME(CopyNodesToElementsNodes)

	const NodeData &from;
	ParameterData &to;

	CopyNodesToElementsNodes(const NodeData &from, ParameterData &to): from(from), to(to)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void now();
};

struct AverageElementsNodesToNodes: public OperatorBuilder {
	GET_NAME(AverageElementsNodesToNodes)

	const ParameterData &from;
	NodeData &to;

	AverageElementsNodesToNodes(const ParameterData &from, NodeData &to): from(from), to(to)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void now();
};

struct CopyNodesToBoundaryNodes: public OperatorBuilder {
	GET_NAME(CopyNodesToBoundaryNodes)

	const NodeData &from;
	BoundaryParameterPack &to;

	CopyNodesToBoundaryNodes(const NodeData &from, BoundaryParameterPack &to): from(from), to(to)
	{
		for (size_t r = 0; r < to.regions.size(); ++r) {
			to.regions[r].isset = true;
		}
	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void now();
};

struct CopyBoundaryRegionsSettingToNodes: public OperatorBuilder {
	GET_NAME(CopyBoundaryRegionsSettingToNodes)

	const std::map<std::string, ECFExpression> &from;
	NodeData &to;

	CopyBoundaryRegionsSettingToNodes(const std::map<std::string, ECFExpression> &from, NodeData &to): from(from), to(to)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void now();
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_ */
