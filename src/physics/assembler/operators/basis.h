
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"

#include "physics/assembler/operator.hpp"
#include "physics/assembler/kernels/heattransfer.kernel.opt.h"

namespace espreso {

struct Basis: public ElementOperatorBuilder {
	GET_NAME(Basis)

	bool build(HeatTransferKernelOpt &kernel) override
	{
		kernel.integration.N.resize();
		kernel.integration.dN.resize();
		kernel.integration.weight.resize();
		{
			int index = 0;
			for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei, ++index) {
				Operator::Link(index).inputs(1).outputs(kernel.integration.N, kernel.integration.dN, kernel.integration.weight);
				double *n = (kernel.integration.N.data->begin() + index)->data();
				double *dn = (kernel.integration.dN.data->begin() + index)->data();
				double *w = (kernel.integration.weight.data->begin() + index)->data();
				esint nsize = (*Mesh::edata[ei->code].N)[0].ncols * (*Mesh::edata[ei->code].N)[0].nrows;
				for (size_t gp = 0; gp < Mesh::edata[ei->code].N->size(); ++gp) {
					memcpy(n + gp * nsize, (*Mesh::edata[ei->code].N)[gp].vals, sizeof(double) * nsize);
					memcpy(dn + Mesh::edata[ei->code].dimension * gp * nsize, (*Mesh::edata[ei->code].dN)[gp].vals, sizeof(double) * Mesh::edata[ei->code].dimension * nsize);
				}
				memcpy(w, Mesh::edata[ei->code].weighFactor->data(), sizeof(double) * nsize);
			}
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				kernel.integration.boundary.N.regions[r].resize();
				kernel.integration.boundary.dN.regions[r].resize();
				kernel.integration.boundary.weight.regions[r].resize();

				int index = 0;
				for (auto ei = info::mesh->boundaryRegions[r]->eintervals.begin(); ei != info::mesh->boundaryRegions[r]->eintervals.end(); ++ei, ++index) {
					Operator::Link(index).inputs(1).outputs(kernel.integration.boundary.N.regions[r], kernel.integration.boundary.dN.regions[r], kernel.integration.boundary.weight.regions[r]);
					double *n = (kernel.integration.boundary.N.regions[r].data->begin() + index)->data();
					double *dn = (kernel.integration.boundary.dN.regions[r].data->begin() + index)->data();
					double *w = (kernel.integration.boundary.weight.regions[r].data->begin() + index)->data();
					esint nsize = (*Mesh::edata[ei->code].N)[0].ncols * (*Mesh::edata[ei->code].N)[0].nrows;
					for (size_t gp = 0; gp < Mesh::edata[ei->code].N->size(); ++gp) {
						memcpy(n + gp * nsize, (*Mesh::edata[ei->code].N)[gp].vals, sizeof(double) * nsize);
						memcpy(dn + Mesh::edata[ei->code].dimension * gp * nsize, (*Mesh::edata[ei->code].dN)[gp].vals, sizeof(double) * Mesh::edata[ei->code].dimension * nsize);
					}
					memcpy(w, Mesh::edata[ei->code].weighFactor->data(), sizeof(double) * nsize);
				}
			}
		}
		return true;
	}

	void apply(int interval)
	{

	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_BASIS_H_ */
