
#include "uniformbuilder.direct.h"

#include "basis/containers/serializededata.h"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/nodestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

UniformBuilderDirectPattern::UniformBuilderDirectPattern(HeatTransferLoadStepConfiguration &configuration)
{
    dofs = 1;
    shape = Matrix_Shape::UPPER;
    type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    for (auto mat = info::ecf->heat_transfer.material_set.begin(); mat != info::ecf->heat_transfer.material_set.end(); ++mat) {
        if (info::ecf->heat_transfer.materials.find(mat->second)->second.thermal_conductivity.model == ThermalConductivityConfiguration::MODEL::ANISOTROPIC) {
            shape = Matrix_Shape::FULL;
            type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
        }
    }
    if (configuration.translation_motions.size()) {
        shape = Matrix_Shape::FULL;
        type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
    }

    std::vector<esint> indices;
    bregion.resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (configuration.temperature.find(region->name) != configuration.temperature.end()) {
            bregion[r].dirichlet = 1;
        }
    }

    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (bregion[r].dirichlet) {
            for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
                for (int d = 0; d < dofs; ++d) {
                    indices.push_back(*n * dofs + d);
                }
            }
        }
    }
    bregion[0].b = indices; // use the first region to store indices permutation;
    utils::sortAndRemoveDuplicates(indices);
    bregion[0].indices = indices;
    for (size_t i = 0; i < bregion[0].b.size(); ++i) {
        bregion[0].b[i] = std::lower_bound(indices.begin(), indices.end(), bregion[0].b[i]) - indices.begin();
    }
    buildPattern(dofs);
}

UniformBuilderDirectPattern::UniformBuilderDirectPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = info::mesh->dimension * multiplicity;
    shape = Matrix_Shape::UPPER;
    type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;

    std::vector<esint> indices;
    bregion.resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        auto expr = configuration.displacement.find(region->name);
        if (expr != configuration.displacement.end()) {
            for (int m = 0; m < multiplicity; ++m) {
                for (int d = 0; d < dofs; ++d) {
                    if (expr->second.data[d].isset) {
                        bregion[r].dirichlet += 1 << (dofs * m + d);
                    }
                }
            }
        }
    }

    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
            for (int m = 0; m < multiplicity; ++m) {
                for (int d = 0; d < dofs; ++d) {
                    if (bregion[r].dirichlet & (1 << (dofs * m + d))) {
                        indices.push_back(*n * dofs * multiplicity + m * dofs + d);
                    }
                }
            }
        }
    }
    bregion[0].b = indices; // use the first region to store indices permutation;
    utils::sortAndRemoveDuplicates(indices);
    bregion[0].indices = indices;
    for (size_t i = 0; i < bregion[0].b.size(); ++i) {
        bregion[0].b[i] = std::lower_bound(indices.begin(), indices.end(), bregion[0].b[i]) - indices.begin();
    }
    buildPattern(dofs * multiplicity);
}

void UniformBuilderDirectPattern::buildPattern(int dofs)
{
    double start = eslog::time();
    eslog::info(" == LINEAR SYSTEM                                                               DISTRIBUTED == \n");
    eslog::info(" == DOFS PER NODE                                                                         %d == \n", dofs);

    std::vector<std::vector<esint> > sSize(info::mesh->neighbors.size()), rSize(info::mesh->neighbors.size());
    std::vector<esint> begin(info::mesh->nodes->size + 1, 1); // add diagonal
    for (auto enodes = info::mesh->elements->nodes->cbegin(); enodes != info::mesh->elements->nodes->cend(); ++enodes) {
        for (auto n = enodes->begin(); n != enodes->end(); ++n) {
            begin[*n] += enodes->size() - 1; // do not count diagonal
        }
    }

    auto ranks = info::mesh->nodes->ranks->cbegin();
    for (esint n = 0; n < info::mesh->nodes->size; ++n, ++ranks) {
        auto neigh = info::mesh->neighbors.begin();
        for (auto r = ranks->begin(); r != ranks->end(); ++r) {
            if (*r != info::mpi::rank) {
                while (*neigh < *r) { ++neigh; }
                sSize[neigh - info::mesh->neighbors.begin()].push_back(begin[n] - 1);
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sSize, rSize, info::mesh->neighbors)) {
        eslog::internalFailure("send size of node intervals.\n");
    }
    sSize.clear();

    ranks = info::mesh->nodes->ranks->cbegin();
    std::vector<esint> rIndex(info::mesh->neighbors.size());
    for (esint n = 0; n < info::mesh->nodes->size; ++n, ++ranks) {
        auto neigh = info::mesh->neighbors.begin();
        for (auto r = ranks->begin(); r != ranks->end(); ++r) {
            if (*r != info::mpi::rank) {
                while (*neigh < *r) { ++neigh; }
                begin[n] += rSize[neigh - info::mesh->neighbors.begin()][rIndex[neigh - info::mesh->neighbors.begin()]++];
            }
        }
    }

    utils::sizesToOffsets(begin);
    std::vector<esint> end = begin;
    std::vector<esint, initless_allocator<esint> > indices(begin.back());
    for (esint n = 0; n < info::mesh->nodes->size; ++n) {
        indices[end[n]++] = info::mesh->nodes->uniqInfo.position[n]; // diagonal
    }

    for (auto enodes = info::mesh->elements->nodes->cbegin(); enodes != info::mesh->elements->nodes->cend(); ++enodes) {
        for (auto from = enodes->begin(); from != enodes->end(); ++from) {
            for (auto to = enodes->begin(); to != enodes->end(); ++to) {
                if (*from != *to) {
                    indices[end[*from]++] = info::mesh->nodes->uniqInfo.position[*to];
                }
            }
        }
    }

    std::vector<std::vector<esint> > sIndices(info::mesh->neighbors.size()), rIndices(info::mesh->neighbors.size());
    ranks = info::mesh->nodes->ranks->cbegin();
    for (esint n = 0; n < info::mesh->nodes->size; ++n, ++ranks) {
        auto neigh = info::mesh->neighbors.begin();
        for (auto r = ranks->begin(); r != ranks->end(); ++r) {
            if (*r != info::mpi::rank) {
                while (*neigh < *r) { ++neigh; }
                for (auto i = begin[n] + 1; i < end[n]; ++i) {
                    sIndices[neigh - info::mesh->neighbors.begin()].push_back(indices[i]);
                }
            }
        }
    }

    if (!Communication::exchangeUnknownSize(sIndices, rIndices, info::mesh->neighbors)) {
        eslog::internalFailure("send size of node intervals.\n");
    }
    sIndices.clear();

    ranks = info::mesh->nodes->ranks->cbegin();
    std::fill(rIndex.begin(), rIndex.end(), 0);
    std::vector<esint> rDataIndex(info::mesh->neighbors.size());
    for (esint n = 0; n < info::mesh->nodes->size; ++n, ++ranks) {
        auto neigh = info::mesh->neighbors.begin();
        for (auto r = ranks->begin(); r != ranks->end(); ++r) {
            if (*r != info::mpi::rank) {
                while (*neigh < *r) { ++neigh; }
                size_t nn = neigh - info::mesh->neighbors.begin();
                for (auto i = 0; i < rSize[nn][rIndex[nn]]; ++i) {
                    indices[end[n]++] = rIndices[nn][rDataIndex[nn]++];
                }
                ++rIndex[nn];
            }
        }
    }

    size_t count = 0;
    #pragma omp parallel for reduction(+:count)
    for (esint n = 0; n < info::mesh->nodes->size; ++n) {
        std::sort(indices.begin() + begin[n], indices.begin() + end[n]);
        esint unique = begin[n];
        for (auto i = begin[n] + 1; i < end[n]; ++i) {
            if (indices[unique] != indices[i]) {
                indices[++unique] = indices[i];
            }
        }
        end[n] = unique + 1;
        count += end[n] - begin[n];
    }
    count *= dofs * dofs;
    eslog::info(" == NON-ZERO VALUES                                                          %14lu == \n", count);

    elements.nrows = dofs * (info::mesh->nodes->uniqInfo.nhalo + info::mesh->nodes->uniqInfo.size);
    elements.ncols = dofs * info::mesh->nodes->uniqInfo.totalSize;
    elements.row.reserve(dofs * count);
    elements.column.reserve(dofs * count);
    elements.A.reserve(dofs * dofs * indices.size());
    elements.b.reserve(dofs * info::mesh->elements->nodes->datatarray().size());

    std::vector<esint, initless_allocator<esint> > offset;
    offset.reserve(info::mesh->nodes->size);
    for (esint n = 0, size = 0; n < info::mesh->nodes->size; ++n) {
        offset.push_back(size);
        for (int r = 0; r < dofs; ++r) {
            for (auto i = begin[n]; i < end[n]; ++i) {
                for (int c = 0; c < dofs; ++c) {
                    elements.row.push_back(info::mesh->nodes->uniqInfo.position[n] * dofs + r);
                    elements.column.push_back(indices[i] * dofs + c);
                }
            }
        }
        size += dofs * dofs * (end[n] - begin[n]);
    }
    for (auto enodes = info::mesh->elements->nodes->cbegin(); enodes != info::mesh->elements->nodes->cend(); ++enodes) {
        for (int rd = 0; rd < dofs; ++rd) {
            for (auto from = enodes->begin(); from != enodes->end(); ++from) {
                auto ibegin = indices.begin() + begin[*from];
                auto iend = indices.begin() + end[*from];
                esint roffset = offset[*from] + rd * dofs * (iend - ibegin);
                for (int cd = 0; cd < dofs; ++cd) {
                    for (auto to = enodes->begin(); to != enodes->end(); ++to) {
                        elements.A.push_back(roffset + dofs * (std::lower_bound(ibegin, iend, info::mesh->nodes->uniqInfo.position[*to]) - ibegin) + cd);
                    }
                }
                elements.b.push_back(*from * dofs + rd);
            }
        }
    }

    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            for (auto e = info::mesh->boundaryRegions[r]->elements->cbegin(); e != info::mesh->boundaryRegions[r]->elements->cend(); ++e) {
                for (int rd = 0; rd < dofs; ++rd) {
                    for (auto from = e->begin(); from != e->end(); ++from) {
                        auto ibegin = indices.begin() + begin[*from];
                        auto iend = indices.begin() + end[*from];
                        esint roffset = offset[*from] + rd * dofs * (iend - ibegin);
                        for (int cd = 0; cd < dofs; ++cd) {
                            for (auto to = e->begin(); to != e->end(); ++to) {
                                bregion[r].A.push_back(roffset + dofs * (std::lower_bound(ibegin, iend, info::mesh->nodes->uniqInfo.position[*to]) - ibegin) + cd);
                            }
                        }
                        bregion[r].b.push_back(*from * dofs + rd);
                    }
                }
            }
        } else {
            bregion[r].b.reserve(info::mesh->boundaryRegions[r]->nodes->datatarray().size());
            for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().cbegin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); ++n) {
                for (int rd = 0; rd < dofs; ++rd) {
                    bregion[r].b.push_back(*n * dofs + rd);
                }
            }
        }
    }

    decomposition.begin = dofs * info::mesh->nodes->uniqInfo.offset;
    decomposition.end = dofs * (info::mesh->nodes->uniqInfo.offset + info::mesh->nodes->uniqInfo.size);
    decomposition.totalSize = dofs * info::mesh->nodes->uniqInfo.totalSize;

    decomposition.neighbors = info::mesh->neighbors;
    decomposition.neighDOF.resize(decomposition.neighbors.size() + 1, decomposition.begin); // the last is my offset
    decomposition.halo.clear();
    decomposition.halo.reserve(dofs * info::mesh->nodes->uniqInfo.nhalo);
    for (esint n = 0; n < info::mesh->nodes->uniqInfo.nhalo; ++n) {
        for (int dof = 0; dof < dofs; ++dof) {
            decomposition.halo.push_back(dofs * info::mesh->nodes->uniqInfo.position[n] + dof);
        }
    }

    std::vector<esint> dBuffer = { decomposition.begin };
    if (!Communication::gatherUniformNeighbors(dBuffer, decomposition.neighDOF, decomposition.neighbors)) {
        eslog::internalFailure("cannot exchange matrix decomposition info.\n");
    }

    size_t nonzeros[2] = { elements.row.size(), bregion[0].b.size() };
    Communication::allReduce(&nonzeros, NULL, 2, MPITools::getType<size_t>().mpitype, MPI_SUM);


    eslog::info(" == DIRICHLET SIZE                                                           %14lu == \n", nonzeros[1]);
    eslog::info(" == LINEAR SYSTEM SIZE                                                       %14d == \n", decomposition.totalSize);
    eslog::info(" == NON-ZERO VALUES                                                          %14lu == \n", nonzeros[0]);
    eslog::info(" == NON-ZERO FILL-IN RATIO                                                         %7.4f%% == \n", 100.0 * nonzeros[0] / decomposition.totalSize / decomposition.totalSize);
    eslog::info(" == COMPOSITION RUNTIME                                                          %8.3f s == \n", eslog::time() - start);
    eslog::info(" ============================================================================================= \n");
}

