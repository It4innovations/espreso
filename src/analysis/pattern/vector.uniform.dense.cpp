
#include "vector.uniform.dense.h"

#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/nodestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

VectorUniformDense::VectorUniformDense(HeatTransferLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = 1;
    buildPattern(dofs);
}

VectorUniformDense::VectorUniformDense(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = info::mesh->dimension * multiplicity;
    buildPattern(dofs);
}

void VectorUniformDense::buildPattern(int dofs)
{
    elements.permutation.reserve(dofs * info::mesh->elements->nodes->datatarray().size());

    for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        elements.offset.push_back(elements.permutation.size());
        const auto &einterval = info::mesh->elements->eintervals[i];
        for (auto enodes = info::mesh->elements->nodes->cbegin() + einterval.begin; enodes != info::mesh->elements->nodes->cbegin() + einterval.end; ++enodes) {
            for (int rd = 0; rd < dofs; ++rd) {
                for (auto from = enodes->begin(); from != enodes->end(); ++from) {
                    elements.permutation.push_back(*from * dofs + rd);
                }
            }
        }
    }

    boundary.resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (info::mesh->boundaryRegions[r]->dimension) {
            for (size_t i = 0; i < region->eintervals.size(); ++i) {
                boundary[r].offset.push_back(boundary[r].permutation.size());
                for (auto e = region->elements->cbegin() + region->eintervals[i].begin; e != region->elements->cbegin() + region->eintervals[i].end; ++e) {
                    for (int rd = 0; rd < dofs; ++rd) {
                        for (auto from = e->begin(); from != e->end(); ++from) {
                            boundary[r].permutation.push_back(*from * dofs + rd);
                        }
                    }
                }
            }
        }
        for (size_t t = 0; t < region->nodes->threads(); ++t) {
            boundary[r].offset.push_back(boundary[r].permutation.size());
            for (auto n = region->nodes->datatarray().cbegin(t); n != region->nodes->datatarray().cend(t); ++n) {
                for (int d = 0; d < dofs; ++d) {
                    boundary[r].permutation.push_back(*n * dofs + d);
                }
            }
        }
    }
}

template <typename T>
void VectorUniformDense::Sync<T>::init(DecompositionDirect &decomposition)
{
    if (info::mpi::size == 1 || neighbors.size()) {
        return;
    }
    neighbors = decomposition.neighbors;
    sBuffer.resize(neighbors.size());
    rBuffer.resize(neighbors.size());
    rOffset.resize(neighbors.size());

    std::vector<std::vector<esint> > sbuffer(neighbors.size()), rbuffer(neighbors.size());

    for (size_t n = 0, r = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
        nOffset.push_back(r);
        while (r < decomposition.halo.size() && decomposition.halo[r] < decomposition.neighDOF[n + 1]) {
            sbuffer[n].push_back(decomposition.halo[r++]);
        }
        sBuffer[n].resize(r - nOffset.back());
    }

    if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, neighbors)) {
        eslog::internalFailure("receive MatrixDenseDistributed pattern.\n");
    }

    for (size_t n = 0; n < neighbors.size(); ++n) {
        for (size_t i = 0; i < rbuffer[n].size(); ++i) {
            rOffset[n].push_back(rbuffer[n][i] - decomposition.begin + decomposition.halo.size());
        }
        rBuffer[n].resize(rOffset[n].size());
    }
}

template <typename T>
void VectorUniformDense::Sync<T>::gatherFromUpper(Vector_Distributed<Vector_Dense, T> &v)
{
    for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
        memcpy(sBuffer[n].data(), v.cluster.vals + nOffset[n], sizeof(double) * sBuffer[n].size());
    }

    if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
        eslog::internalFailure("receive MatrixCSRDistribution data.\n");
    }

    for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
        for (size_t i = 0; i < rOffset[n].size(); ++i) {
            v.cluster.vals[rOffset[n][i]] += rBuffer[n][i];
        }
    }
}

template <typename T>
void VectorUniformDense::Sync<T>::scatterToUpper(Vector_Distributed<Vector_Dense, T> &v)
{
    for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
        for (size_t i = 0; i < rOffset[n].size(); ++i) {
            rBuffer[n][i] = v.cluster.vals[rOffset[n][i]];
        }
    }

    if (!Communication::receiveLowerKnownSize(rBuffer, sBuffer, neighbors)) {
        eslog::internalFailure("scatter VectorDenseDistributed data.\n");
    }

    for (size_t n = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
        memcpy(v.cluster.vals + nOffset[n], sBuffer[n].data(), sizeof(double) * sBuffer[n].size());
    }
}

template struct VectorUniformDense::Sync<double>;


