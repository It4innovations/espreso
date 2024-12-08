
#include "matrix.uniform.direct.h"

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

#include <set>

using namespace espreso;

MatrixUniformDirect::MatrixUniformDirect(int DOFs)
{
    dofs = DOFs;
    shape = Matrix_Shape::UPPER;
    type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
    buildPattern(dofs);
}

MatrixUniformDirect::MatrixUniformDirect(HeatTransferLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = multiplicity;
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

    buildPattern(dofs);
}

MatrixUniformDirect::MatrixUniformDirect(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity)
{
    dofs = info::mesh->dimension * multiplicity;
    shape = Matrix_Shape::UPPER;
    type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;

    if (configuration.mode == StructuralMechanicsLoadStepConfiguration::MODE::NONLINEAR) {
        type = Matrix_Type::REAL_SYMMETRIC_INDEFINITE;
    }

//    if (configuration.type == LoadStepSolverConfiguration::TYPE::TRANSIENT) {
//        type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
//    }

    buildPattern(dofs);
}

void MatrixUniformDirect::buildPattern(int dofs)
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

    pattern.nrows = dofs * (info::mesh->nodes->uniqInfo.nhalo + info::mesh->nodes->uniqInfo.size);
    pattern.ncols = dofs * info::mesh->nodes->uniqInfo.totalSize;
    pattern.row.reserve(pattern.nrows + 1);
    pattern.column.reserve(dofs * count);
    elements.permutation.reserve(dofs * dofs * indices.size());

    std::vector<esint, initless_allocator<esint> > offset;
    offset.reserve(info::mesh->nodes->size);
    pattern.row.push_back(0);
    for (esint n = 0, size = 0; n < info::mesh->nodes->size; ++n) {
        offset.push_back(size);
        for (int r = 0; r < dofs; ++r) {
            for (auto i = begin[n]; i < end[n]; ++i) {
                for (int c = 0; c < dofs; ++c) {
                    pattern.column.push_back(indices[i] * dofs + c);
                }
            }
            pattern.row.push_back(pattern.column.size());
        }
        size += dofs * dofs * (end[n] - begin[n]);
    }
    for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
        elements.offset.push_back(elements.permutation.size());
        for (auto enodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin; enodes != info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].end; ++enodes) {
            for (int rd = 0; rd < dofs; ++rd) {
                for (auto from = enodes->begin(); from != enodes->end(); ++from) {
                    auto ibegin = indices.begin() + begin[*from];
                    auto iend = indices.begin() + end[*from];
                    esint roffset = offset[*from] + rd * dofs * (iend - ibegin);
                    for (int cd = 0; cd < dofs; ++cd) {
                        for (auto to = enodes->begin(); to != enodes->end(); ++to) {
                            elements.permutation.push_back(roffset + dofs * (std::lower_bound(ibegin, iend, info::mesh->nodes->uniqInfo.position[*to]) - ibegin) + cd);
                        }
                    }
                }
            }
        }
    }

    boundary.resize(info::mesh->boundaryRegions.size());
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
        if (region->dimension) {
            for (size_t i = 0; i < region->eintervals.size(); ++i) {
                boundary[r].offset.push_back(boundary[r].permutation.size());
                for (auto e = region->elements->cbegin() + region->eintervals[i].begin; e != region->elements->cbegin() + region->eintervals[i].end; ++e) {
                    for (int rd = 0; rd < dofs; ++rd) {
                        for (auto from = e->begin(); from != e->end(); ++from) {
                            auto ibegin = indices.begin() + begin[*from];
                            auto iend = indices.begin() + end[*from];
                            esint roffset = offset[*from] + rd * dofs * (iend - ibegin);
                            for (int cd = 0; cd < dofs; ++cd) {
                                for (auto to = e->begin(); to != e->end(); ++to) {
                                    boundary[r].permutation.push_back(roffset + dofs * (std::lower_bound(ibegin, iend, info::mesh->nodes->uniqInfo.position[*to]) - ibegin) + cd);
                                }
                            }
                        }
                    }
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

    size_t nonzeros = pattern.column.size();
    Communication::allReduce(&nonzeros, NULL, 1, MPITools::getType(nonzeros).mpitype, MPI_SUM);


    eslog::info(" == LINEAR SYSTEM SIZE                                                       %14d == \n", decomposition.totalSize);
    eslog::info(" == NON-ZERO VALUES                                                          %14lu == \n", nonzeros);
    eslog::info(" == NON-ZERO FILL-IN RATIO                                                         %7.4f%% == \n", 100.0 * nonzeros / decomposition.totalSize / decomposition.totalSize);
    eslog::info(" == COMPOSITION RUNTIME                                                          %8.3f s == \n", eslog::time() - start);
    eslog::info(" ============================================================================================= \n");
}

template <typename T>
void MatrixUniformDirect::Sync<T>::init(MatrixUniformDirect &m)
{
    if (info::mpi::size == 1 || neighbors.size()) {
        return;
    }

    neighbors = m.decomposition.neighbors;
    sBuffer.resize(neighbors.size());
    rBuffer.resize(neighbors.size());
    rOffset.resize(neighbors.size());

    std::vector<std::vector<esint> > sbuffer(neighbors.size()), rbuffer(neighbors.size());

    auto halo = m.decomposition.halo.begin();
    for (size_t n = 0, r = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
        size_t size = 0;
        nOffset.push_back(m.pattern.row[r] - Indexing::CSR);
        while (halo != m.decomposition.halo.end() && *halo < m.decomposition.neighDOF[n + 1]) {
            sbuffer[n].push_back(*halo);
            sbuffer[n].push_back(m.pattern.row[r + 1] - m.pattern.row[r]);
            size += m.pattern.row[r + 1] - m.pattern.row[r];
            for (esint c = m.pattern.row[r]; c < m.pattern.row[r + 1]; ++c) {
                sbuffer[n].push_back(m.pattern.column[c - Indexing::CSR]);
            }
            ++halo; ++r;
        }
        sBuffer[n].resize(size);
    }

    if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, neighbors)) {
        eslog::internalFailure("receive MatrixCSRDistribution pattern.\n");
    }

    for (size_t n = 0; n < neighbors.size(); ++n) {
        for (size_t i = 0; i < rbuffer[n].size(); ) {
            esint *c = m.pattern.column.data() + m.pattern.row[rbuffer[n][i++] - m.decomposition.begin + m.decomposition.halo.size()] - Indexing::CSR;
            esint columns = rbuffer[n][i++];
            for (esint cc = 0; cc < columns; ++cc, ++i) {
                while (*c != rbuffer[n][i]) { ++c; }
                rOffset[n].push_back(c - m.pattern.column.data());
            }
        }
        rBuffer[n].resize(rOffset[n].size());
    }
}

template <typename T>
void MatrixUniformDirect::Sync<T>::gatherFromUpper(Matrix_Distributed<T> &m)
{
    for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
        memcpy(sBuffer[n].data(), m.cluster.vals + nOffset[n], sizeof(double) * sBuffer[n].size());
    }

    if (!Communication::receiveUpperUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
        eslog::internalFailure("receive MatrixCSRDistribution data.\n");
    }

    for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
        for (size_t i = 0; i < rOffset[n].size(); ++i) {
            m.cluster.vals[rOffset[n][i]] += rBuffer[n][i];
        }
    }
}

template <typename T>
void MatrixUniformDirect::Sync<T>::scatterToUpper(Matrix_Distributed<T> &m)
{
    for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
        for (size_t i = 0; i < rOffset[n].size(); ++i) {
            rBuffer[n][i] = m.cluster.vals[rOffset[n][i]];
        }
    }

    if (!Communication::receiveUpperUnknownSize(rBuffer, sBuffer, info::mesh->neighbors)) {
        eslog::internalFailure("receive MatrixCSRDistribution data.\n");
    }

    for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
        memcpy(m.cluster.vals + nOffset[n], sBuffer[n].data(), sizeof(double) * sBuffer[n].size());
    }
}

template <typename T>
void MatrixUniformDirect::Apply<T>::init(MatrixUniformDirect &m)
{
    esint nhalo = m.decomposition.halo.size();
    esint msize = m.decomposition.end - m.decomposition.begin;

    { // build neighbors of neighbors
        std::vector<std::pair<esint, esint> > nn;
        std::vector<std::vector<std::pair<esint, esint> > > rnn(m.decomposition.neighbors.size());
        for (size_t n = 0; n < m.decomposition.neighbors.size(); ++n) {
            nn.push_back(std::make_pair(m.decomposition.neighbors[n], m.decomposition.neighDOF[n]));
        }

        if (!Communication::exchangeUnknownSize(nn, rnn, m.decomposition.neighbors)) {
            eslog::internalFailure("cannot exchange neighbors neighbors.\n");
        }
        size_t upper = 0;
        for (size_t n = 0; n < rnn.size(); ++n) {
            upper += rnn[n].size();
        }
        nn.reserve(nn.size() + upper);
        for (size_t n = 0; n < rnn.size(); ++n) {
            nn.insert(nn.end(), rnn[n].begin(), rnn[n].end());
        }
        utils::sortAndRemoveDuplicates(nn);

        for (size_t n = 0; n < nn.size(); ++n) {
            if (nn[n].first != info::mpi::rank) {
                neighbors.push_back(nn[n].first);
                nDOF.push_back(nn[n].second);
            }
        }
    }

    sBuffer.resize(neighbors.size());
    rBuffer.resize(neighbors.size());
    sOffset.resize(neighbors.size());
    rOffset.resize(neighbors.size());

    std::set<esint> lower, higher;
    std::vector<std::set<esint> > send(neighbors.size());
    for (esint r = nhalo; r < m.pattern.nrows; ++r) {
        for (esint c = m.pattern.row[r] - Indexing::CSR; c < m.pattern.row[r + 1] - Indexing::CSR; ++c) {
            if (m.pattern.column[c] - Indexing::CSR < m.decomposition.begin) {
                lower.insert(m.pattern.column[c] - Indexing::CSR);
                auto n = std::lower_bound(nDOF.begin(), nDOF.end(), m.pattern.column[c] - Indexing::CSR + 1) - nDOF.begin() - 1;
                send[n].insert(r);
            }
            if (m.decomposition.end <= m.pattern.column[c] - Indexing::CSR) {
                higher.insert(m.pattern.column[c] - Indexing::CSR);
                auto n = std::lower_bound(nDOF.begin(), nDOF.end(), m.pattern.column[c] - Indexing::CSR + 1) - nDOF.begin() - 1;
                send[n].insert(r);
            }
        }
    }

    size_t ni = 0, ci = 0;
    for (auto c = lower.begin(); c != lower.end(); ++c, ++ci) {
        while (ni + 1 < nDOF.size() && nDOF[ni + 1] <= *c) {
            ++ni;
        }
        rOffset[ni].push_back(ci);
    }
    for (auto c = higher.begin(); c != higher.end(); ++c, ++ci) {
        while (ni + 1 < nDOF.size() && nDOF[ni + 1] <= *c) {
            ++ni;
        }
        rOffset[ni].push_back(ci + msize);
    }
    for (size_t n = 0; n < rOffset.size(); ++n) {
        rBuffer[n].resize(rOffset[n].size());
        sBuffer[n].resize(send[n].size());
        sOffset[n].insert(sOffset[n].end(), send[n].begin(), send[n].end());
    }

    localM.type = Matrix_Type::COMPLEX_STRUCTURALLY_SYMMETRIC;
    localM.shape = m.shape;
    localM._allocated.nrows = localM.nrows = msize;
    localM._allocated.ncols = localM.ncols = msize + ci;
    localM._allocated.nnz   = localM.nnz   = m.pattern.row[m.pattern.nrows] - m.pattern.row[nhalo];
    localM._allocated.rows  = localM.rows  = localM.ator.template allocate<esint>(localM.nrows + 1);
    localM._allocated.cols  = localM.cols  = localM.ator.template allocate<esint>(localM.nnz);

    offset = lower.size();
    localV.resize(localM.ncols);

    localM.rows[0] = Indexing::CSR;
    for (esint r = nhalo, i = 1; r < m.pattern.nrows; ++r, ++i) {
        localM.rows[i] = localM.rows[i - 1] + m.pattern.row[r + 1] - m.pattern.row[r];
    }

    std::vector<esint> _lower(lower.begin(), lower.end()), _higher(higher.begin(), higher.end());
    for (esint i = 0, c = m.pattern.row[nhalo] - Indexing::CSR; i < localM.nnz; ++i, ++c) {
        if (m.pattern.column[c] - Indexing::CSR < m.decomposition.begin) {
            localM.cols[i] = std::lower_bound(_lower.begin(), _lower.end(), m.pattern.column[c] - Indexing::CSR) - _lower.begin() + Indexing::CSR;
        } else if (m.decomposition.end <= m.pattern.column[c] - Indexing::CSR) {
            localM.cols[i] = std::lower_bound(_higher.begin(), _higher.end(), m.pattern.column[c] - Indexing::CSR) - _higher.begin() + _lower.size() + msize + Indexing::CSR;
        } else {
            localM.cols[i] = m.pattern.column[c] - m.decomposition.begin + _lower.size();
        }
    }
}

template <typename T>
void MatrixUniformDirect::Apply<T>::apply(Matrix_Distributed<T> &m, Vector_Distributed<Vector_Dense, T> &y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> &x)
{
    localM.vals = m.cluster.vals + m.cluster.rows[m.decomposition->halo.size()] - Indexing::CSR;
    spblas.insert(localM);

    for (size_t n = 0; n < sOffset.size(); ++n) {
        for (size_t i = 0; i < sOffset[n].size(); ++i) {
            sBuffer[n][i] = x.cluster.vals[sOffset[n][i]];
        }
    }
    if (!Communication::exchangeKnownSize(sBuffer, rBuffer, neighbors)) {
        eslog::internalFailure("Cannot exchange apply data.\n");
    }
    for (size_t n = 0; n < rOffset.size(); ++n) {
        for (size_t i = 0; i < rOffset[n].size(); ++i) {
            localV.vals[rOffset[n][i]] = rBuffer[n][i];
        }
    }
    std::copy(x.cluster.vals + x.decomposition->halo.size(), x.cluster.vals + x.cluster.size, localV.vals + offset);

    Vector_Dense<T, esint> v;
    v.size = localM.nrows;
    v.vals = y.cluster.vals + m.decomposition->halo.size();
    spblas.apply(v, alpha, beta, localV);
    y.scatter();
}

template struct MatrixUniformDirect::Sync<double>;
template struct MatrixUniformDirect::Apply<double>;



