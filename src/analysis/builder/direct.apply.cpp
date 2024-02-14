
#include "direct.apply.h"
#include "analysis/math/matrix_distributed.h"
#include "analysis/math/vector_distributed.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "wrappers/mpi/communication.h"

#include <set>
#include <vector>

namespace espreso {

template <typename T>
void _buildNeighNeighs(Matrix_CSR_Apply<T> *data, Matrix_Distributed<T> &m)
{
    std::vector<std::pair<esint, esint> > nn;
    std::vector<std::vector<std::pair<esint, esint> > > rnn(m.decomposition->neighbors.size());
    for (size_t n = 0; n < m.decomposition->neighbors.size(); ++n) {
        nn.push_back(std::make_pair(m.decomposition->neighbors[n], m.decomposition->neighDOF[n]));
    }

    if (!Communication::exchangeUnknownSize(nn, rnn, m.decomposition->neighbors)) {
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
            data->neighbors.push_back(nn[n].first);
            data->nDOF.push_back(nn[n].second);
        }
    }
}

template <typename T>
void Matrix_CSR_Apply<T>::init(Matrix_Distributed<T> &m)
{
    const Matrix_CSR<T, esint> &matrix = m.cluster;
    esint nhalo = m.decomposition->halo.size();
    esint msize = m.decomposition->end - m.decomposition->begin;

    _buildNeighNeighs(this, m);

    sBuffer.resize(neighbors.size());
    rBuffer.resize(neighbors.size());
    sOffset.resize(neighbors.size());
    rOffset.resize(neighbors.size());

    std::set<esint> lower, higher;
    std::vector<std::set<esint> > send(neighbors.size());
    for (esint r = nhalo; r < matrix.nrows; ++r) {
        for (esint c = matrix.rows[r] - Indexing::CSR; c < matrix.rows[r + 1] - Indexing::CSR; ++c) {
            if (matrix.cols[c] - Indexing::CSR < m.decomposition->begin) {
                lower.insert(matrix.cols[c] - Indexing::CSR);
                auto n = std::lower_bound(nDOF.begin(), nDOF.end(), matrix.cols[c] - Indexing::CSR + 1) - nDOF.begin() - 1;
                send[n].insert(r);
            }
            if (m.decomposition->end <= matrix.cols[c] - Indexing::CSR) {
                higher.insert(matrix.cols[c] - Indexing::CSR);
                auto n = std::lower_bound(nDOF.begin(), nDOF.end(), matrix.cols[c] - Indexing::CSR + 1) - nDOF.begin() - 1;
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

    localM.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
    localM.shape = m.cluster.shape;
    localM._allocated.nrows = localM.nrows = msize;
    localM._allocated.ncols = localM.ncols = msize + ci;
    localM._allocated.nnz   = localM.nnz   = matrix.rows[matrix.nrows] - matrix.rows[nhalo];
    localM._allocated.rows  = localM.rows  = new esint[localM.nrows + 1];
    localM._allocated.cols  = localM.cols  = new esint[localM.nnz];

    offset = lower.size();
    localV.resize(localM.ncols);

    localM.rows[0] = Indexing::CSR;
    for (esint r = nhalo, i = 1; r < matrix.nrows; ++r, ++i) {
        localM.rows[i] = localM.rows[i - 1] + matrix.rows[r + 1] - matrix.rows[r];
    }

    std::vector<esint> _lower(lower.begin(), lower.end()), _higher(higher.begin(), higher.end());
    for (esint i = 0, c = matrix.rows[nhalo] - Indexing::CSR; i < localM.nnz; ++i, ++c) {
        if (matrix.cols[c] - Indexing::CSR < m.decomposition->begin) {
            localM.cols[i] = std::lower_bound(_lower.begin(), _lower.end(), matrix.cols[c] - Indexing::CSR) - _lower.begin() + Indexing::CSR;
        } else if (m.decomposition->end <= matrix.cols[c] - Indexing::CSR) {
            localM.cols[i] = std::lower_bound(_higher.begin(), _higher.end(), matrix.cols[c] - Indexing::CSR) - _higher.begin() + _lower.size() + msize + Indexing::CSR;
        } else {
            localM.cols[i] = matrix.cols[c] - m.decomposition->begin + _lower.size();
        }
    }
}

template <typename T>
void Matrix_CSR_Apply<T>::apply(Matrix_Distributed<T> &m, Vector_Distributed<Vector_Dense, T> *y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> *x)
{
    localM.vals = m.cluster.vals + m.cluster.rows[m.decomposition->halo.size()] - Indexing::CSR;
    spblas.insert(localM);

    for (size_t n = 0; n < sOffset.size(); ++n) {
        for (size_t i = 0; i < sOffset[n].size(); ++i) {
            sBuffer[n][i] = x->cluster.vals[sOffset[n][i]];
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
    std::copy(x->cluster.vals + x->decomposition->halo.size(), x->cluster.vals + x->cluster.size, localV.vals + offset);

    Vector_Dense<T, esint> v;
    v.size = localM.nrows;
    v.vals = y->cluster.vals + m.decomposition->halo.size();
    spblas.apply(v, alpha, beta, localV);
    y->scatter();
}

template struct Matrix_CSR_Apply<double>;

}
