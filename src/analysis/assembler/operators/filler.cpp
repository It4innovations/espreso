
#include "filler.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/heattransfer.generator.h"

#include "analysis/scheme/steadystate.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

using namespace espreso;

void HeatTransfer::connect(SteadyState &scheme)
{
	switch (scheme.K->shape) {
	case Matrix_Shape::FULL:  generateElementOperators<GeneralMatricFiller>(etype, elementOps, 1, elements.stiffness, scheme.K); break;
	case Matrix_Shape::UPPER: generateElementOperators<SymmetricMatricFiller>(etype, elementOps, 1, elements.stiffness, scheme.K); break;
	}

	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		switch (info::mesh->dimension) {
		case 2:
			switch (info::mesh->boundaryRegions[r]->dimension) {
			case 0:
				for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
					if (configuration.temperature.end() == configuration.temperature.find(info::mesh->boundaryRegions[r]->name)) {
						// sector filler
					} else {
						if (boundaryOps[r][t].size()) {
							boundaryOps[r][t].push_back(generateNodeSetter2D<VectorSetter>(r, t, 1, scheme.dirichlet,
									[] (auto &element, const size_t &n, const size_t &dof, const size_t &s) -> double& { return element.temperature[n][s]; }));
						}
					}
				}
				break;
			case 1:
				break;
			}
			break;
		case 3:
			switch (info::mesh->boundaryRegions[r]->dimension) {
			case 0:
				for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
					if (configuration.temperature.end() == configuration.temperature.find(info::mesh->boundaryRegions[r]->name)) {
						// sector filler
					} else {
						if (boundaryOps[r][t].size()) {
							boundaryOps[r][t].push_back(generateNodeSetter3D<VectorSetter>(r, t, 1, scheme.dirichlet,
									[] (auto &element, const size_t &n, const size_t &dof, const size_t &s) -> double& { return element.temperature[n][s]; }));
						}
					}
				}
				break;
			case 1:
				break;
			case 2:
				break;
			}
			break;
		}
	}
}
