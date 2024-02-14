
#include "direct.synchronization.h"
#include "analysis/math/matrix_distributed.h"
#include "analysis/math/vector_distributed.h"
#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_csr.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "wrappers/mpi/communication.h"

#include <cstring>

namespace espreso {

template <typename T>
void Matrix_CSR_Sync<T>::init(Matrix_Distributed<T> &m)
{
    if (info::mpi::size == 1 || neighbors.size()) {
        return;
    }

    neighbors = m.decomposition->neighbors;
    sBuffer.resize(neighbors.size());
    rBuffer.resize(neighbors.size());
    rOffset.resize(neighbors.size());

    std::vector<std::vector<esint> > sbuffer(neighbors.size()), rbuffer(neighbors.size());

    const Matrix_CSR<T> &csr = m.cluster;

    auto halo = m.decomposition->halo.begin();
    for (size_t n = 0, r = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
        size_t size = 0;
        nOffset.push_back(csr.rows[r] - Indexing::CSR);
        while (halo != m.decomposition->halo.end() && *halo < m.decomposition->neighDOF[n + 1]) {
            sbuffer[n].push_back(*halo);
            sbuffer[n].push_back(csr.rows[r + 1] - csr.rows[r]);
            size += csr.rows[r + 1] - csr.rows[r];
            for (esint c = csr.rows[r]; c < csr.rows[r + 1]; ++c) {
                sbuffer[n].push_back(csr.cols[c - Indexing::CSR]);
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
            esint *c = csr.cols + csr.rows[rbuffer[n][i++] - m.decomposition->begin + m.decomposition->halo.size()] - Indexing::CSR;
            esint columns = rbuffer[n][i++];
            for (esint cc = 0; cc < columns; ++cc, ++i) {
                while (*c != rbuffer[n][i]) { ++c; }
                rOffset[n].push_back(c - csr.cols);
            }
        }
        rBuffer[n].resize(rOffset[n].size());
    }
}

template <typename T>
void Matrix_CSR_Sync<T>::gatherFromUpper(Matrix_Distributed<T> &m)
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
void Matrix_CSR_Sync<T>::scatterToUpper(Matrix_Distributed<T> &m)
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
void Vector_Dense_Sync<T>::init(Vector_Distributed<Vector_Dense, T> &v)
{
    if (info::mpi::size == 1 || neighbors.size()) {
        return;
    }
    neighbors = v.decomposition->neighbors;
    sBuffer.resize(neighbors.size());
    rBuffer.resize(neighbors.size());
    rOffset.resize(neighbors.size());

    std::vector<std::vector<esint> > sbuffer(neighbors.size()), rbuffer(neighbors.size());

    for (size_t n = 0, r = 0; n < neighbors.size() && neighbors[n] < info::mpi::rank; ++n) {
        nOffset.push_back(r);
        while (r < v.decomposition->halo.size() && v.decomposition->halo[r] < v.decomposition->neighDOF[n + 1]) {
            sbuffer[n].push_back(v.decomposition->halo[r++]);
        }
        sBuffer[n].resize(r - nOffset.back());
    }

    if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, neighbors)) {
        eslog::internalFailure("receive MatrixDenseDistributed pattern.\n");
    }

    for (size_t n = 0; n < neighbors.size(); ++n) {
        for (size_t i = 0; i < rbuffer[n].size(); ++i) {
            rOffset[n].push_back(rbuffer[n][i] - v.decomposition->begin + v.decomposition->halo.size());
        }
        rBuffer[n].resize(rOffset[n].size());
    }
}

template <typename T>
void Vector_Dense_Sync<T>::gatherFromUpper(Vector_Distributed<Vector_Dense, T> &v)
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
void Vector_Dense_Sync<T>::scatterToUpper(Vector_Distributed<Vector_Dense, T> &v)
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

template <typename T>
void Vector_Sparse_Sync<T>::init(Vector_Distributed<Vector_Sparse, T> &v)
{

}

template <typename T>
void Vector_Sparse_Sync<T>::gatherFromUpper(Vector_Distributed<Vector_Sparse, T> &v)
{

}

template <typename T>
void Vector_Sparse_Sync<T>::scatterToUpper(Vector_Distributed<Vector_Sparse, T> &v)
{

}

template struct Matrix_CSR_Sync<double>;
template struct Vector_Dense_Sync<double>;
template struct Vector_Sparse_Sync<double>;

}


