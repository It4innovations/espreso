
#include "synchronization.h"

#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "math/primitives/matrix_csr.h"
#include "math/generalization/matrix_distributed.h"
#include "wrappers/mpi/communication.h"

#include <cstring>
#include <complex>

namespace espreso {

template <typename T>
void _init(Data_Synchronization<Matrix_CSR, T> *sync, Matrix_Distributed<Matrix_CSR, T> &m)
{
	if (info::mpi::size == 1 || sync->neighbors.size()) {
		return;
	}

	sync->neighbors = m.distribution->neighbors;
	sync->sBuffer.resize(sync->neighbors.size());
	sync->rBuffer.resize(sync->neighbors.size());
	sync->rOffset.resize(sync->neighbors.size());

	std::vector<std::vector<esint> > sbuffer(sync->neighbors.size()), rbuffer(sync->neighbors.size());

	const Matrix_CSR<T> &csr = m.cluster;

	auto halo = m.distribution->halo.begin();
	for (size_t n = 0, r = 0; n < sync->neighbors.size() && sync->neighbors[n] < info::mpi::rank; ++n) {
		size_t size = 0;
		sync->nOffset.push_back(csr.rows[r] - Indexing::CSR);
		while (halo != m.distribution->halo.end() && *halo < m.distribution->neighDOF[n + 1]) {
			sbuffer[n].push_back(*halo);
			sbuffer[n].push_back(csr.rows[r + 1] - csr.rows[r]);
			size += csr.rows[r + 1] - csr.rows[r];
			for (esint c = csr.rows[r]; c < csr.rows[r + 1]; ++c) {
				sbuffer[n].push_back(csr.cols[c - Indexing::CSR]);
			}
			++halo; ++r;
		}
		sync->sBuffer[n].resize(size);
	}

	if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, sync->neighbors)) {
		eslog::internalFailure("receive MatrixCSRDistribution pattern.\n");
	}

	for (size_t n = 0; n < sync->neighbors.size(); ++n) {
		for (size_t i = 0; i < rbuffer[n].size(); ) {
			esint *c = csr.cols + csr.rows[rbuffer[n][i++] - m.distribution->begin + m.distribution->halo.size()] - Indexing::CSR;
			esint columns = rbuffer[n][i++];
			for (esint cc = 0; cc < columns; ++cc, ++i) {
				while (*c != rbuffer[n][i]) { ++c; }
				sync->rOffset[n].push_back(c - csr.cols);
			}
		}
		sync->rBuffer[n].resize(sync->rOffset[n].size());
	}
}

template <typename T>
void _gatherFromUpper(Data_Synchronization<Matrix_CSR, T> *sync, Matrix_Distributed<Matrix_CSR, T> &m)
{
	for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
		memcpy(sync->sBuffer[n].data(), m.cluster.vals + sync->nOffset[n], sizeof(double) * sync->sBuffer[n].size());
	}

	if (!Communication::receiveUpperUnknownSize(sync->sBuffer, sync->rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("receive MatrixCSRDistribution data.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		for (size_t i = 0; i < sync->rOffset[n].size(); ++i) {
			m.cluster.vals[sync->rOffset[n][i]] += sync->rBuffer[n][i];
		}
	}
}

template <typename T>
void _scatterToUpper(Data_Synchronization<Matrix_CSR, T> *sync, Matrix_Distributed<Matrix_CSR, T> &m)
{
	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		for (size_t i = 0; i < sync->rOffset[n].size(); ++i) {
			sync->rBuffer[n][i] = m.cluster.vals[sync->rOffset[n][i]];
		}
	}

	if (!Communication::receiveUpperUnknownSize(sync->rBuffer, sync->sBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("receive MatrixCSRDistribution data.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
		memcpy(m.cluster.vals + sync->nOffset[n], sync->sBuffer[n].data(), sizeof(double) * sync->sBuffer[n].size());
	}
}

template<> void Data_Synchronization<Matrix_CSR, double>::init(Matrix_Distributed<Matrix_CSR, double> &m) { _init(this, m); }
template<> void Data_Synchronization<Matrix_CSR, std::complex<double> >::init(Matrix_Distributed<Matrix_CSR, std::complex<double> > &m) { _init(this, m); }

template<> void Data_Synchronization<Matrix_CSR, double>::gatherFromUpper(Matrix_Distributed<Matrix_CSR, double> &m) { _gatherFromUpper(this, m); }
template<> void Data_Synchronization<Matrix_CSR, std::complex<double> >::gatherFromUpper(Matrix_Distributed<Matrix_CSR, std::complex<double> > &m) { _gatherFromUpper(this, m); }

template<> void Data_Synchronization<Matrix_CSR, double>::scatterToUpper(Matrix_Distributed<Matrix_CSR, double> &m) { _scatterToUpper(this, m); }
template<> void Data_Synchronization<Matrix_CSR, std::complex<double> >::scatterToUpper(Matrix_Distributed<Matrix_CSR, std::complex<double> > &m) { _scatterToUpper(this, m); }


template <typename T>
void _init(Data_Synchronization<Vector_Dense, T> *sync, Vector_Distributed<Vector_Dense, T> &v)
{
	if (info::mpi::size == 1 || sync->neighbors.size()) {
		return;
	}
	sync->neighbors = v.distribution->neighbors;
	sync->sBuffer.resize(sync->neighbors.size());
	sync->rBuffer.resize(sync->neighbors.size());
	sync->rOffset.resize(sync->neighbors.size());

	std::vector<std::vector<esint> > sbuffer(sync->neighbors.size()), rbuffer(sync->neighbors.size());

	for (size_t n = 0, r = 0; n < sync->neighbors.size() && sync->neighbors[n] < info::mpi::rank; ++n) {
		sync->nOffset.push_back(r);
		while (r < v.distribution->halo.size() && v.distribution->halo[r] < v.distribution->neighDOF[n + 1]) {
			sbuffer[n].push_back(v.distribution->halo[r++]);
		}
		sync->sBuffer[n].resize(r - sync->nOffset.back());
	}

	if (!Communication::receiveUpperUnknownSize(sbuffer, rbuffer, sync->neighbors)) {
		eslog::internalFailure("receive MatrixDenseDistributed pattern.\n");
	}

	for (size_t n = 0; n < sync->neighbors.size(); ++n) {
		for (size_t i = 0; i < rbuffer[n].size(); ++i) {
			sync->rOffset[n].push_back(rbuffer[n][i] - v.distribution->begin + v.distribution->halo.size());
		}
		sync->rBuffer[n].resize(sync->rOffset[n].size());
	}
}

template <typename T>
void _gatherFromUpper(Data_Synchronization<Vector_Dense, T> *sync, Vector_Distributed<Vector_Dense, T> &v)
{
	for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
		memcpy(sync->sBuffer[n].data(), v.cluster.vals + sync->nOffset[n], sizeof(double) * sync->sBuffer[n].size());
	}

	if (!Communication::receiveUpperUnknownSize(sync->sBuffer, sync->rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("receive MatrixCSRDistribution data.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		for (size_t i = 0; i < sync->rOffset[n].size(); ++i) {
			v.cluster.vals[sync->rOffset[n][i]] += sync->rBuffer[n][i];
		}
	}
}

template <typename T>
void _scatterToUpper(Data_Synchronization<Vector_Dense, T> *sync, Vector_Distributed<Vector_Dense, T> &v)
{
	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		for (size_t i = 0; i < sync->rOffset[n].size(); ++i) {
			sync->rBuffer[n][i] = v.cluster.vals[sync->rOffset[n][i]];
		}
	}

	if (!Communication::receiveLowerKnownSize(sync->rBuffer, sync->sBuffer, sync->neighbors)) {
		eslog::internalFailure("scatter VectorDenseDistributed data.\n");
	}

	for (size_t n = 0; n < sync->neighbors.size() && sync->neighbors[n] < info::mpi::rank; ++n) {
		memcpy(v.cluster.vals + sync->nOffset[n], sync->sBuffer[n].data(), sizeof(double) * sync->sBuffer[n].size());
	}
}

template<> void Data_Synchronization<Vector_Dense, double>::init(Vector_Distributed<Vector_Dense, double> &v) { _init(this, v); }
template<> void Data_Synchronization<Vector_Dense, std::complex<double> >::init(Vector_Distributed<Vector_Dense, std::complex<double> > &v) { _init(this, v); }

template<> void Data_Synchronization<Vector_Dense, double>::gatherFromUpper(Vector_Distributed<Vector_Dense, double> &v) { _gatherFromUpper(this, v); }
template<> void Data_Synchronization<Vector_Dense, std::complex<double> >::gatherFromUpper(Vector_Distributed<Vector_Dense, std::complex<double> > &v) { _gatherFromUpper(this, v); }

template<> void Data_Synchronization<Vector_Dense, double>::scatterToUpper(Vector_Distributed<Vector_Dense, double> &v) { _scatterToUpper(this, v); }
template<> void Data_Synchronization<Vector_Dense, std::complex<double> >::scatterToUpper(Vector_Distributed<Vector_Dense, std::complex<double> > &v) { _scatterToUpper(this, v); }

}


