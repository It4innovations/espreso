
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "mesh/store/nodestore.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t dimension>
struct CopyParameter: public ActionOperator {
	CopyParameter(int interval, const ParameterData &from, ParameterData &to)
	: from(from, interval),
	  to(to, interval)
	{

	}

	InputParameterIterator from;
	OutputParameterIterator to;

	void operator++()
	{
		++from;
		++to;
	}

	void move(int n)
	{
		from += n;
		to += n;
	}

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			for (size_t d = 0; d < dimension; ++d) {
				to.data[gpindex * dimension + d] = from.data[gpindex * dimension + d];
			}
		}
	}
};

template <size_t dimension>
struct CopyNodesToEnodes: public ActionOperator {
	CopyNodesToEnodes(int interval, const NodeData &from, serializededata<esint, esint>::const_iterator procNodes, ParameterData &to)
	: from(from),
	  procNodes(procNodes),
	  to(to, interval)
	{

	}

	const NodeData &from;
	serializededata<esint, esint>::const_iterator procNodes;
	OutputParameterIterator to;

	void operator++()
	{
		++procNodes;
		++to;
	}

	void move(int n)
	{
		procNodes += n;
		to += n;
	}

	void operator()()
	{
		size_t i = 0;
		for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++i) {
			for (size_t d = 0; d < dimension; ++d) {
				to[i * dimension + d] = from.data[*n * dimension + d];
			}
		}
	}
};

////////////////

//
//class NodeData;
//class ECFExpression;
//
//template <class TParent>
//struct CopyParameters: public TParent {
//	const ParameterData &from;
//	ParameterData &to;
//
//	CopyParameters(AX_HeatTransfer &module, const ParameterData &from, ParameterData &to, const char *name)
//	: TParent(name), assembler(assembler), from(from), to(to)
//	{
//		to.addInput(from);
//		to.resize();
//		module.addParameter(to);
//	}
//
//	void apply(int interval)
//	{
//		if (to.update[interval]) {
//			if (from.isconst[interval]) {
//				auto indata = from.data->begin() + interval;
//				auto outdata = to.data->begin() + interval;
//				for (size_t i = 0; i < outdata->size() / indata->size(); ++i) {
//					memcpy(outdata->data() + i * indata->size(), indata->data(), sizeof(double) * indata->size());
//				}
//			} else {
//				memcpy((to.data->begin() + interval)->data(), (from.data->begin() + interval)->data(), sizeof(double) * (from.data->begin() + interval)->size());
//			}
//		}
//	}
//};
//
//struct CopyElementParameters: public CopyParameters<ElementOperatorBuilder> {
//	using CopyParameters<ElementOperatorBuilder>::CopyParameters;
//};
//
//struct CopyBoundaryParameters: public CopyParameters<BoundaryOperatorBuilder> {
//	using CopyParameters<BoundaryOperatorBuilder>::CopyParameters;
//};
//
//struct CopyNodesToBoundaryNodes: public OperatorBuilder {
//	const NodeData &from;
//	BoundaryParameterPack &to;
//
//	CopyNodesToBoundaryNodes(AX_HeatTransfer &module, const NodeData &from, BoundaryParameterPack &to, const char *name): OperatorBuilder(name), from(from), to(to)
//	{
//		for (size_t r = 0; r < to.regions.size(); ++r) {
//			to.regions[r].isset = true;
//		}
//		to.addInput(&from);
//		to.resize();
//		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
//			module.addParameter(to.regions[r]);
//		}
//	}
//
//	void now()
//	{
//		if (Operator::print > 1) printf("%s\n", name);
//		#pragma omp parallel for
//		for (int t = 0; t < info::env::threads; t++) {
//			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
//				if (info::mesh->boundaryRegions[r]->dimension) {
//					for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; ++d) {
//						for (esint ii = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; ii < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++ii) {
//							auto i = (to.regions[r].data->begin() + ii)->data();
//							auto procNodes = info::mesh->boundaryRegions[r]->nodes->begin() + info::mesh->boundaryRegions[r]->eintervals[ii].begin;
//							for (auto e = info::mesh->boundaryRegions[r]->eintervals[ii].begin; e != info::mesh->boundaryRegions[r]->eintervals[ii].end; ++e, ++procNodes) {
//								for (auto n = procNodes->begin(); n != procNodes->end(); ++n) {
//									for (int dim = 0; dim < from.dimension; ++dim, ++i) {
//										*i = from.data[*n * from.dimension + dim];
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//};
//
//struct CopyBoundaryRegionsSettingToNodes: public OperatorBuilder {
//	const std::map<std::string, ECFExpression> &from;
//	NodeData &to;
//
//	CopyBoundaryRegionsSettingToNodes(const std::map<std::string, ECFExpression> &from, NodeData &to, const char *name): OperatorBuilder(name), from(from), to(to)
//	{
//
//	}
//
//	void now()
//	{
//		if (Operator::print > 1) printf("%s\n", name);
//		#pragma omp parallel for
//		for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
//			for (auto it = from.begin(); it != from.end(); ++it) {
//				serializededata<esint, esint>* nodes = info::mesh->bregion(it->first)->nodes;
//				it->second.evaluator->evalSelectedDense(nodes->datatarray().size(t), nodes->datatarray().begin(t), it->second.evaluator->params, to.data.data());
//			}
//		}
//	}
//};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_ */
