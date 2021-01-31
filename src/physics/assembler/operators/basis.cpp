
#include "basis.h"

#include "esinfo/meshinfo.h"
#include "math/matrix.dense.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

using namespace espreso;

bool Basis::build(HeatTransferModuleOpt &kernel)
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
			esint nodes = Mesh::edata[ei->code].nodes;
			esint gps = Mesh::edata[ei->code].weighFactor->size();
			for (esint gp = 0; gp < gps; ++gp) {
				memcpy(n + gp * nodes, (*Mesh::edata[ei->code].N)[gp].vals, sizeof(double) * nodes);
				memcpy(dn + Mesh::edata[ei->code].dimension * gp * nodes, (*Mesh::edata[ei->code].dN)[gp].vals, sizeof(double) * Mesh::edata[ei->code].dimension * nodes);
			}
			memcpy(w, Mesh::edata[ei->code].weighFactor->data(), sizeof(double) * gps);
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
				esint nodes = Mesh::edata[ei->code].nodes;
				esint gps = Mesh::edata[ei->code].weighFactor->size();
				for (esint gp = 0; gp < gps; ++gp) {
					memcpy(n + gp * nodes, (*Mesh::edata[ei->code].N)[gp].vals, sizeof(double) * nodes);
					memcpy(dn + Mesh::edata[ei->code].dimension * gp * nodes, (*Mesh::edata[ei->code].dN)[gp].vals, sizeof(double) * Mesh::edata[ei->code].dimension * nodes);
				}
				memcpy(w, Mesh::edata[ei->code].weighFactor->data(), sizeof(double) * gps);
			}
		}
	}
	return true;
}
