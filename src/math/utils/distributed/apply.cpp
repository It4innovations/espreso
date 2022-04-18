
#include "apply.h"

#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "math/math.h"
#include "math/primitives/matrix_csr.h"
#include "math/physics/matrix_distributed.h"
#include "wrappers/mpi/communication.h"

#include <complex>
#include <set>
#include <algorithm>

namespace espreso {

template <typename T>
void _buildNeighNeighs(Data_Apply<Matrix_CSR, T> *data, Matrix_Distributed<Matrix_CSR, T> &m)
{
	std::vector<std::pair<esint, esint> > nn;
	std::vector<std::vector<std::pair<esint, esint> > > rnn(m.distribution->neighbors.size());
	for (size_t n = 0; n < m.distribution->neighbors.size(); ++n) {
		nn.push_back(std::make_pair(m.distribution->neighbors[n], m.distribution->neighDOF[n]));
	}

	if (!Communication::exchangeUnknownSize(nn, rnn, m.distribution->neighbors)) {
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
void _init(Data_Apply<Matrix_CSR, T> *data, Matrix_Distributed<Matrix_CSR, T> &m)
{
	const Matrix_CSR<T> &matrix = m.cluster;
	esint nhalo = m.distribution->halo.size();
	esint msize = m.distribution->end - m.distribution->begin;

	_buildNeighNeighs(data, m);

	data->sBuffer.resize(data->neighbors.size());
	data->rBuffer.resize(data->neighbors.size());
	data->sOffset.resize(data->neighbors.size());
	data->rOffset.resize(data->neighbors.size());

	std::set<esint> lower, higher;
	std::vector<std::set<esint> > send(data->neighbors.size());
	for (esint r = nhalo; r < matrix.nrows; ++r) {
		for (esint c = matrix.rows[r] - Indexing::CSR; c < matrix.rows[r + 1] - Indexing::CSR; ++c) {
			if (matrix.cols[c] - Indexing::CSR < m.distribution->begin) {
				lower.insert(matrix.cols[c] - Indexing::CSR);
				auto n = std::lower_bound(data->nDOF.begin(), data->nDOF.end(), matrix.cols[c] - Indexing::CSR + 1) - data->nDOF.begin() - 1;
				send[n].insert(r);
			}
			if (m.distribution->end <= matrix.cols[c] - Indexing::CSR) {
				higher.insert(matrix.cols[c] - Indexing::CSR);
				auto n = std::lower_bound(data->nDOF.begin(), data->nDOF.end(), matrix.cols[c] - Indexing::CSR + 1) - data->nDOF.begin() - 1;
				send[n].insert(r);
			}
		}
	}

	size_t ni = 0, ci = 0;
	for (auto c = lower.begin(); c != lower.end(); ++c, ++ci) {
		while (ni + 1 < data->nDOF.size() && data->nDOF[ni + 1] <= *c) {
			++ni;
		}
		data->rOffset[ni].push_back(ci);
	}
	for (auto c = higher.begin(); c != higher.end(); ++c, ++ci) {
		while (ni + 1 < data->nDOF.size() && data->nDOF[ni + 1] <= *c) {
			++ni;
		}
		data->rOffset[ni].push_back(ci + msize);
	}
	for (size_t n = 0; n < data->rOffset.size(); ++n) {
		data->rBuffer[n].resize(data->rOffset[n].size());
		data->sBuffer[n].resize(send[n].size());
		data->sOffset[n].insert(data->sOffset[n].end(), send[n].begin(), send[n].end());
	}

	data->m.type = Matrix_Type::REAL_STRUCTURALLY_SYMMETRIC;
	data->m.shape = m.cluster.shape;
	data->m._allocated.nrows = data->m.nrows = msize;
	data->m._allocated.ncols = data->m.ncols = msize + ci;
	data->m._allocated.nnz = data->m.nnz = matrix.rows[matrix.nrows] - matrix.rows[nhalo];
	data->m._allocated.rows = data->m.rows = new esint[data->m.nrows + 1];
	data->m._allocated.cols = data->m.cols = new esint[data->m.nnz];

	data->offset = lower.size();
	data->v.size = data->m.ncols;
	data->v._allocated.vals = data->v.vals = new T[data->m.ncols];

	data->m.rows[0] = Indexing::CSR;
	for (esint r = nhalo, i = 1; r < matrix.nrows; ++r, ++i) {
		data->m.rows[i] = data->m.rows[i - 1] + matrix.rows[r + 1] - matrix.rows[r];
	}

	std::vector<esint> _lower(lower.begin(), lower.end()), _higher(higher.begin(), higher.end());
	for (esint i = 0, c = matrix.rows[nhalo] - Indexing::CSR; i < data->m.nnz; ++i, ++c) {
		if (matrix.cols[c] - Indexing::CSR < m.distribution->begin) {
			data->m.cols[i] = std::lower_bound(_lower.begin(), _lower.end(), matrix.cols[c] - Indexing::CSR) - _lower.begin() + Indexing::CSR;
		} else if (m.distribution->end <= matrix.cols[c] - Indexing::CSR) {
			data->m.cols[i] = std::lower_bound(_higher.begin(), _higher.end(), matrix.cols[c] - Indexing::CSR) - _higher.begin() + _lower.size() + msize + Indexing::CSR;
		} else {
			data->m.cols[i] = matrix.cols[c] - m.distribution->begin + _lower.size();
		}
	}
}

template <typename T>
void _commit(Data_Apply<Matrix_CSR, T> *data, Matrix_Distributed<Matrix_CSR, T> &m)
{
	data->m.vals = m.cluster.vals + m.cluster.rows[m.distribution->halo.size()] - Indexing::CSR;
	math::commit(data->m);
}

template <typename T>
void _apply(Data_Apply<Matrix_CSR, T> *data, Vector_Distributed<Vector_Dense, T> *y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> *x)
{
	for (size_t n = 0; n < data->sOffset.size(); ++n) {
		for (size_t i = 0; i < data->sOffset[n].size(); ++i) {
			data->sBuffer[n][i] = x->cluster.vals[data->sOffset[n][i]];
		}
	}
	if (!Communication::exchangeKnownSize(data->sBuffer, data->rBuffer, data->neighbors)) {
		eslog::internalFailure("Cannot exchange apply data.\n");
	}
	for (size_t n = 0; n < data->rOffset.size(); ++n) {
		for (size_t i = 0; i < data->rOffset[n].size(); ++i) {
			data->v.vals[data->rOffset[n][i]] = data->rBuffer[n][i];
		}
	}
	std::copy(x->cluster.vals + x->distribution->halo.size(), x->cluster.vals + x->cluster.size, data->v.vals + data->offset);

	Vector_Dense<T> v;
	v.size = data->m.nrows;
	v.vals = y->cluster.vals + y->distribution->halo.size();
	math::apply(v, alpha, data->m, beta, data->v);
	y->scatter();
}

template<> void Data_Apply<Matrix_CSR, double>::init(Matrix_Distributed<Matrix_CSR, double> &m) { _init(this, m); }
template<> void Data_Apply<Matrix_CSR, std::complex<double> >::init(Matrix_Distributed<Matrix_CSR, std::complex<double> > &m) { _init(this, m); }

template<> void Data_Apply<Matrix_CSR, double>::commit(Matrix_Distributed<Matrix_CSR, double> &m) { _commit(this, m); }
template<> void Data_Apply<Matrix_CSR, std::complex<double> >::commit(Matrix_Distributed<Matrix_CSR, std::complex<double> > &m) { _commit(this, m); }

template<> void Data_Apply<Matrix_CSR, double>::apply(Vector_Distributed<Vector_Dense, double> *y, const double &alpha, const double &beta, const Vector_Distributed<Vector_Dense, double> *x) { _apply(this, y, alpha, beta, x); }
template<> void Data_Apply<Matrix_CSR, std::complex<double> >::apply(Vector_Distributed<Vector_Dense, std::complex<double>> *y, const std::complex<double> &alpha, const std::complex<double> &beta, const Vector_Distributed<Vector_Dense, std::complex<double> > *x) { _apply(this, y, alpha, beta, x); }

}
