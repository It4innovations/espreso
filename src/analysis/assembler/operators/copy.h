
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "basis/evaluator/evaluator.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "wrappers/mpi/communication.h"

#include <cstdio>
#include <map>

namespace espreso {

class NodeData;
class ECFExpression;

template <class Assembler, class TParent>
struct CopyParameters: public TParent {
	Assembler &assembler;
	const ParameterData &from;
	ParameterData &to;

	CopyParameters(Assembler &assembler, const ParameterData &from, ParameterData &to, const char *name)
	: TParent(name), assembler(assembler), from(from), to(to)
	{

	}

	bool build() override
	{
		to.addInput(from);
		to.resize();
		assembler.addParameter(to);
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

template <class Assembler>
struct CopyElementParameters: public CopyParameters<Assembler, ElementOperatorBuilder> {
	using CopyParameters<Assembler, ElementOperatorBuilder>::CopyParameters;
};

template <class Assembler>
struct CopyBoundaryParameters: public CopyParameters<Assembler, BoundaryOperatorBuilder> {
	using CopyParameters<Assembler, BoundaryOperatorBuilder>::CopyParameters;
};

template <typename Assembler>
struct CopyNodesToElementsNodes: public OperatorBuilder {
	Assembler &assembler;
	const NodeData &from;
	ParameterData &to;

	CopyNodesToElementsNodes(Assembler &assembler, const NodeData &from, ParameterData &to, const char *name): OperatorBuilder(name), assembler(assembler), from(from), to(to)
	{

	}

	bool build() override
	{
		to.addInput(&from);
		to.resize();
		assembler.addParameter(to);
		return true;
	}

	void now()
	{
		if (Operator::print > 1) printf("%s\n", name);
		#pragma omp parallel for
		for (int t = 0; t < info::env::threads; t++) {
			for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
				for (esint ii = info::mesh->elements->eintervalsDistribution[d]; ii < info::mesh->elements->eintervalsDistribution[d + 1]; ++ii) {
					auto i = (to.data->begin() + ii)->data();
					auto procNodes = info::mesh->elements->nodes->begin() + info::mesh->elements->eintervals[ii].begin;
					for (esint e = info::mesh->elements->eintervals[ii].begin; e < info::mesh->elements->eintervals[ii].end; ++e, ++procNodes) {
						for (auto n = procNodes->begin(); n != procNodes->end(); ++n) {
							for (int dim = 0; dim < from.dimension; ++dim, ++i) {
								*i = from.data[*n * from.dimension + dim];
							}
						}
					}
				}
			}
		}
	}
};

template <typename Assembler>
struct AverageElementsNodesToNodes: public OperatorBuilder {
	Assembler &assembler;
	const ParameterData &from;
	NodeData &to;

	AverageElementsNodesToNodes(Assembler &assembler, const ParameterData &from, NodeData &to, const char *name): OperatorBuilder(name), assembler(assembler), from(from), to(to)
	{

	}

	bool build() override
	{
		return true;
	}

	void now()
	{
		if (Operator::print > 1) printf("%s\n", name);
		std::fill(to.data.begin(), to.data.end(), 0);
		auto procNodes = info::mesh->elements->nodes->begin();
		for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			InputParameterIterator input(from, i, to.dimension);
			for (esint e = info::mesh->elements->eintervals[i].begin; e < info::mesh->elements->eintervals[i].end; ++e, ++procNodes) {
				for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++input) {
					for (int d = 0; d < to.dimension; ++d) {
						to.data[*n * to.dimension + d] += input.data[d];
					}
				}
			}
		}

		auto &nelements = info::mesh->nodes->elements->boundarytarray();
		for (size_t i = 0; i < to.data.size() / to.dimension; i++) {
			for (int d = 0; d < to.dimension; d++) {
				to.data[i * to.dimension + d] /= nelements[i + 1] - nelements[i];
			}
		}

		std::vector<std::vector<double> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());

		auto nranks = info::mesh->nodes->ranks->begin();
		for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
			if (nranks->size() > 1) {
				esint noffset = 0;
				for (auto r = nranks->begin(); r != nranks->end(); ++r) {
					if (*r != info::mpi::rank) {
						while (info::mesh->neighbors[noffset] < *r) {
							++noffset;
						}
						for (int d = 0; d < to.dimension; d++) {
							sBuffer[noffset].push_back(to.data[n * to.dimension + d]);
						}
					}
				}
			}
		}

		for (size_t n = 0; n < info::mesh->neighbors.size(); n++) {
			rBuffer[n].resize(sBuffer[n].size());
		}

		if (!Communication::exchangeKnownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
			eslog::internalFailure("exchange diagonal values.\n");
		}

		nranks = info::mesh->nodes->ranks->begin();
		std::vector<esint> nindex(info::mesh->neighbors.size());
		for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
			if (nranks->size() > 1) {
				esint noffset = 0;
				for (auto r = nranks->begin(); r != nranks->end(); ++r) {
					if (*r != info::mpi::rank) {
						while (info::mesh->neighbors[noffset] < *r) {
							++noffset;
						}
						for (int d = 0; d < to.dimension; d++) {
							to.data[n * to.dimension + d] += rBuffer[noffset][nindex[noffset]++];
						}
					}
				}
			}
		}
	}
};

template <typename Assembler>
struct CopyNodesToBoundaryNodes: public OperatorBuilder {
	Assembler &assembler;
	const NodeData &from;
	BoundaryParameterPack &to;

	CopyNodesToBoundaryNodes(Assembler &assembler, const NodeData &from, BoundaryParameterPack &to, const char *name): OperatorBuilder(name), assembler(assembler), from(from), to(to)
	{
		for (size_t r = 0; r < to.regions.size(); ++r) {
			to.regions[r].isset = true;
		}
	}

	bool build() override
	{
		to.addInput(&from);
		to.resize();
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			assembler.addParameter(to.regions[r]);
		}
		return true;
	}

	void now()
	{
		if (Operator::print > 1) printf("%s\n", name);
		#pragma omp parallel for
		for (int t = 0; t < info::env::threads; t++) {
			for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; ++d) {
						for (esint ii = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; ii < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++ii) {
							auto i = (to.regions[r].data->begin() + ii)->data();
							auto procNodes = info::mesh->boundaryRegions[r]->nodes->begin() + info::mesh->boundaryRegions[r]->eintervals[ii].begin;
							for (auto e = info::mesh->boundaryRegions[r]->eintervals[ii].begin; e != info::mesh->boundaryRegions[r]->eintervals[ii].end; ++e, ++procNodes) {
								for (auto n = procNodes->begin(); n != procNodes->end(); ++n) {
									for (int dim = 0; dim < from.dimension; ++dim, ++i) {
										*i = from.data[*n * from.dimension + dim];
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

template <typename Assembler>
struct CopyBoundaryRegionsSettingToNodes: public OperatorBuilder {
	Assembler &assembler;
	const std::map<std::string, ECFExpression> &from;
	NodeData &to;

	CopyBoundaryRegionsSettingToNodes(Assembler &assembler, const std::map<std::string, ECFExpression> &from, NodeData &to, const char *name): OperatorBuilder(name), assembler(assembler), from(from), to(to)
	{

	}

	bool build() override
	{
		return true;
	}

	void now()
	{
		if (Operator::print > 1) printf("%s\n", name);
		#pragma omp parallel for
		for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
			for (auto it = from.begin(); it != from.end(); ++it) {
				serializededata<esint, esint>* nodes = info::mesh->bregion(it->first)->nodes;
				it->second.evaluator->evalSelectedDense(nodes->datatarray().size(t), nodes->datatarray().begin(t), it->second.evaluator->params, to.data.data());
			}
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COPY_H_ */
