
#include "vector.uniform.dense.h"

#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/nodestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

VectorUniformDense::VectorUniformDense(int DOFs)
{
    dofs = DOFs;
    buildPattern(dofs);
}

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
void VectorUniformDense::SyncData<T>::init(DecompositionDirect &decomposition)
{
    if (info::mpi::size == 1 || this->neighbors.size()) {
        return;
    }
    this->neighbors = decomposition.neighbors;
    this->sBuffer.resize(this->neighbors.size());
    this->rBuffer.resize(this->neighbors.size());
    this->rOffset.resize(this->neighbors.size());

    std::vector<std::vector<esint> > sbuffer(this->neighbors.size()), rbuffer(this->neighbors.size());

    for (size_t n = 0, r = 0; n < this->neighbors.size() && this->neighbors[n] < info::mpi::rank; ++n) {
        this->nOffset.push_back(r);
        while (r < decomposition.halo.size() && decomposition.halo[r] < decomposition.neighDOF[n + 1]) {
            sbuffer[n].push_back(decomposition.halo[r++]);
        }
        this->sBuffer[n].resize(r - this->nOffset.back());
    }

    if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, this->neighbors)) {
        eslog::internalFailure("receive MatrixDenseDistributed pattern.\n");
    }

    for (size_t n = 0; n < this->neighbors.size(); ++n) {
        for (size_t i = 0; i < rbuffer[n].size(); ++i) {
            this->rOffset[n].push_back(rbuffer[n][i] - decomposition.begin + decomposition.halo.size());
        }
        this->rBuffer[n].resize(this->rOffset[n].size());
    }
}

template <typename T>
void VectorUniformDense::SyncVector<T>::gatherFromUpper(Vector_Distributed<Vector_Dense, T> &v)
{
    for (size_t n = 0; n < this->neighbors.size() && this->neighbors[n] < info::mpi::rank; ++n) {
        memcpy(this->sBuffer[n].data(), v.cluster.vals + this->nOffset[n], sizeof(double) * this->sBuffer[n].size());
    }

    if (!Communication::receiveUpperUnknownSize(this->sBuffer, this->rBuffer, this->neighbors)) {
        eslog::internalFailure("receive MatrixCSRDistribution data.\n");
    }

    for (size_t n = 0; n < this->neighbors.size(); ++n) {
        for (size_t i = 0; i < this->rOffset[n].size(); ++i) {
            v.cluster.vals[this->rOffset[n][i]] += this->rBuffer[n][i];
        }
    }
}

template <typename T>
void VectorUniformDense::SyncVector<T>::scatterToUpper(Vector_Distributed<Vector_Dense, T> &v)
{
    for (size_t n = 0; n < this->neighbors.size(); ++n) {
        for (size_t i = 0; i < this->rOffset[n].size(); ++i) {
            this->rBuffer[n][i] = v.cluster.vals[this->rOffset[n][i]];
        }
    }

    if (!Communication::receiveLowerKnownSize(this->rBuffer, this->sBuffer, this->neighbors)) {
        eslog::internalFailure("scatter VectorDenseDistributed data.\n");
    }

    for (size_t n = 0; n < this->neighbors.size() && this->neighbors[n] < info::mpi::rank; ++n) {
        memcpy(v.cluster.vals + this->nOffset[n], this->sBuffer[n].data(), sizeof(double) * this->sBuffer[n].size());
    }
}

template <typename T>
void VectorUniformDense::SyncMatrix<T>::gatherFromUpper(Vector_Distributed<Matrix_Dense, T> &v)
{
    for (size_t n = 0; n < this->neighbors.size() && this->neighbors[n] < info::mpi::rank; ++n) {
        memcpy(this->sBuffer[n].data(), v.cluster.vals + this->nOffset[n], sizeof(double) * this->sBuffer[n].size());
    }

    if (!Communication::receiveUpperUnknownSize(this->sBuffer, this->rBuffer, this->neighbors)) {
        eslog::internalFailure("receive MatrixCSRDistribution data.\n");
    }

    for (size_t n = 0; n < this->neighbors.size(); ++n) {
        for (size_t i = 0; i < this->rOffset[n].size(); ++i) {
            v.cluster.vals[this->rOffset[n][i]] += this->rBuffer[n][i];
        }
    }
}

template <typename T>
void VectorUniformDense::SyncMatrix<T>::scatterToUpper(Vector_Distributed<Matrix_Dense, T> &v)
{
    for (size_t n = 0; n < this->neighbors.size(); ++n) {
        for (size_t i = 0; i < this->rOffset[n].size(); ++i) {
            this->rBuffer[n][i] = v.cluster.vals[this->rOffset[n][i]];
        }
    }

    if (!Communication::receiveLowerKnownSize(this->rBuffer, this->sBuffer, this->neighbors)) {
        eslog::internalFailure("scatter VectorDenseDistributed data.\n");
    }

    for (size_t n = 0; n < this->neighbors.size() && this->neighbors[n] < info::mpi::rank; ++n) {
        memcpy(v.cluster.vals + this->nOffset[n], this->sBuffer[n].data(), sizeof(double) * this->sBuffer[n].size());
    }
}

template struct VectorUniformDense::SyncData<double>;
template struct VectorUniformDense::SyncVector<double>;
template struct VectorUniformDense::SyncMatrix<double>;



