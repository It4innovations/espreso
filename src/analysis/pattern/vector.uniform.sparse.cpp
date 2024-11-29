
#include "vector.uniform.sparse.h"

#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"

#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/nodestore.h"

using namespace espreso;

VectorUniformSparse::VectorUniformSparse(int DOFs)
{
    dirichlet.resize(info::mesh->boundaryRegions.size());
    std::vector<int> isset(DOFs * info::mesh->nodes->size);
    buildPattern(isset);
}

VectorUniformSparse::VectorUniformSparse(HeatTransferLoadStepConfiguration &configuration, int multiplicity)
{
    int dofs = multiplicity;
    dirichlet.resize(info::mesh->boundaryRegions.size());
    std::vector<int> isset(dofs * info::mesh->nodes->size);
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (configuration.temperature.find(region->name) != configuration.temperature.end()) {
            dirichlet[r].filter = 1;
            for (size_t t = 0; t < region->nodes->threads(); ++t) {
                dirichlet[r].offset.push_back(dirichlet[r].permutation.size());
                for (auto n = region->nodes->datatarray().cbegin(t); n != region->nodes->datatarray().cend(t); ++n) {
                    for (int d = 0; d < dofs; ++d) {
                        dirichlet[r].permutation.push_back(*n * dofs + d);
                        isset[*n * dofs + d] = 1;
                    }
                }
            }
        }
    }
    buildPattern(isset);
}

VectorUniformSparse::VectorUniformSparse(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity)
{
    int dofs = info::mesh->dimension * multiplicity;
    dirichlet.resize(info::mesh->boundaryRegions.size());
    std::vector<int> isset(dofs * info::mesh->nodes->size);
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        auto expr = configuration.displacement.find(region->name);
        if (expr != configuration.displacement.end()) {
            for (int m = 0; m < multiplicity; ++m) {
                for (int d = 0; d < dofs; ++d) {
                    if (expr->second.data[d].isset) {
                        dirichlet[r].filter += 1 << (dofs * m + d);
                    }
                }
            }
            for (size_t t = 0; t < region->nodes->threads(); ++t) {
                dirichlet[r].offset.push_back(dirichlet[r].permutation.size());
                for (auto n = region->nodes->datatarray().cbegin(t); n != region->nodes->datatarray().cend(t); ++n) {
                    for (int m = 0; m < multiplicity; ++m) {
                        for (int d = 0; d < dofs; ++d) {
                            if (expr->second.data[d].isset) {
                                if (dirichlet[r].filter & (1 << (dofs * m + d))) {
                                    dirichlet[r].permutation.push_back(*n * dofs * multiplicity + m * dofs + d);
                                    isset[*n * dofs * multiplicity + m * dofs + d] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    buildPattern(isset);
}

void VectorUniformSparse::buildPattern(std::vector<int> &isset)
{
    for (size_t i = 0, offset = 0; i < isset.size(); ++i) {
        if (isset[i]) {
            indices.push_back(i);
            isset[i] = offset++;
        }
    }
    for (size_t r = 0; r < dirichlet.size(); ++r) {
        for (size_t i = 0; i < dirichlet[r].permutation.size(); ++i) {
            dirichlet[r].permutation[i] = isset[dirichlet[r].permutation[i]];
        }
    }
}



