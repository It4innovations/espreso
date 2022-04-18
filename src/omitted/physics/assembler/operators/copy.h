
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/boundaryregionstore.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/modules/heattransfer.module.opt.h"

#include <cstdio>
#include <map>

namespace espreso {

class NodeData;
class ECFExpression;

template <class TParent>
struct CopyParameters: public TParent {
	const ParameterData &from;
	ParameterData &to;

	CopyParameters(const ParameterData &from, ParameterData &to, const char *name)
	: TParent(name), from(from), to(to)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		to.addInput(from);
		to.resize();
		kernel.addParameter(to);
		return true;
	}

	void apply(int interval)
	{
		if (to.update[interval]) {
			if (from.isconst[interval]) {
				auto indata = from.data->begin() + interval;
				auto outdata = to.data->begin() + interval;
				for (size_t i = 0; i < outdata->size() / indata->size(); ++i) {
					memcpy(outdata->data() + i * indata->size(), indata->data(), sizeof(double) * indata->size());
				}
			} else {
				memcpy((to.data->begin() + interval)->data(), (from.data->begin() + interval)->data(), sizeof(double) * (from.data->begin() + interval)->size());
			}
		}
	}
};

struct CopyElementParameters: public CopyParameters<ElementOperatorBuilder> {
	using CopyParameters<ElementOperatorBuilder>::CopyParameters;
};
struct CopyBoundaryParameters: public CopyParameters<BoundaryOperatorBuilder> {
	using CopyParameters<BoundaryOperatorBuilder>::CopyParameters;
};

struct CopyNodesToElementsNodes: public OperatorBuilder {
	const NodeData &from;
	ParameterData &to;

	CopyNodesToElementsNodes(const NodeData &from, ParameterData &to, const char *name): OperatorBuilder(name), from(from), to(to)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		to.addInput(&from);
		to.resize();
		kernel.addParameter(to);
		return true;
	}

	void now();
};

struct AverageElementsNodesToNodes: public OperatorBuilder {
	const ParameterData &from;
	NodeData &to;

	AverageElementsNodesToNodes(const ParameterData &from, NodeData &to, const char *name): OperatorBuilder(name), from(from), to(to)
	{

	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		return true;
	}

	void now();
};

struct CopyNodesToBoundaryNodes: public OperatorBuilder {
	const NodeData &from;
	BoundaryParameterPack &to;

	CopyNodesToBoundaryNodes(const NodeData &from, BoundaryParameterPack &to, const char *name): OperatorBuilder(name), from(from), to(to)
	{
		for (size_t r = 0; r < to.regions.size(); ++r) {
			to.regions[r].isset = true;
		}
	}

	bool build(HeatTransferModuleOpt &kernel) override
	{
		to.addInput(&from);
		to.resize();
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			kernel.addParameter(to.regions[r]);
		}
		return true;
	}

	void now();
};

struct CopyBoundaryRegionsSettingToNodes: public OperatorBuilder {
	const std::map<std::string, ECFExpression> &from;
	NodeData &to;

	CopyBoundaryRegionsSettingToNodes(const std::map<std::string, ECFExpression> &from, NodeData &to, const char *name): OperatorBuilder(name), from(from), to(to)
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
