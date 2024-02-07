
#include "regularization.elasticity.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "math/wrappers/math.blas.h"
#include "math/wrappers/math.lapack.h"

namespace espreso {

template struct RegularizationElasticity<double>;

template <typename T>
static void getNtNNtN(Matrix_Dense<T> &N, Matrix_Dense<T> &NtNNtN)
{
	Matrix_Dense<T> _N(N), NNt;
	math::blas::AAt(N, NNt);
	math::lapack::solve(NNt, _N);
	math::blas::multiply(T{1}, N, _N, T{0}, NtNNtN, true);
}

template <typename TFix, typename T>
static void setRegMat(std::vector<TFix> &fix, Matrix_CSR<T> &RegMat, esint size)
{
	RegMat.resize(size, size, (fix.size() - 1) * fix.size() / 2 + fix.size());
	RegMat.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
	RegMat.shape = Matrix_Shape::UPPER;

	RegMat.rows[0] = Indexing::CSR;
	esint r = 0;
	for (size_t i = 0; i < fix.size(); ++i, ++r) {
		while (r < fix[i].col) {
			RegMat.rows[r + 1] = RegMat.rows[r];
			++r;
		}
		RegMat.rows[r + 1] = RegMat.rows[r] + fix.size() - i;
		for (size_t j = i; j < fix.size(); ++j) {
			RegMat.cols[RegMat.rows[r] + j - i - Indexing::CSR] = fix[j].col + Indexing::CSR;
		}
	}
	while (r < RegMat.nrows) {
		RegMat.rows[r + 1] = RegMat.rows[r];
		++r;
	}
}

// RegMat = Nt * ((N * Nt)^-1 * N)
//
//    per Fix-Point
// N = [  1  0
//        0  1
//       -y  x ]
//

template <typename T>
void RegularizationElasticity<T>::set2D(esint domain)
{
	std::vector<esint> fixPoints;
	this->getFixPoints(fixPoints, domain);

	const Matrix_CSR<T> &K = feti.K[domain];
	Matrix_Dense<T> &R = feti.R1[domain];
	Matrix_CSR<T> &RegMat = feti.RegMat[domain];
	const FETIDecomposition *decomposition = feti.decomposition;

	struct __fix__ { int col; double value[3]; }; // reorder DOFs to get sorted output
	std::vector<__fix__> fix(2 * fixPoints.size());
	for (size_t i = 0; i < fixPoints.size(); ++i) {
		auto dmap0 = decomposition->dmap->begin() + (2 * fixPoints[i] + 0);
		auto dmap1 = decomposition->dmap->begin() + (2 * fixPoints[i] + 1);
		for (auto di = dmap0->begin(); di != dmap0->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				fix[2 * i + 0].col = di->index;
				fix[2 * i + 0].value[0] = 1;
				fix[2 * i + 0].value[1] = 0;
				fix[2 * i + 0].value[2] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].y;
			}
		}
		for (auto di = dmap1->begin(); di != dmap1->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				fix[2 * i + 1].col = di->index;
				fix[2 * i + 1].value[0] = 0;
				fix[2 * i + 1].value[1] = 1;
				fix[2 * i + 1].value[2] = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].x;
			}
		}
	}
	std::sort(fix.begin(), fix.end(), [] (const __fix__ &f1, const __fix__ &f2) { return f1.col < f2.col; });
	setRegMat(fix, RegMat, K.nrows);

	Matrix_Dense<T> N;
	N.resize(3, fix.size());
	for (size_t i = 0; i < fix.size(); ++i) {
		N.vals[0 * N.ncols + i] = fix[i].value[0];
		N.vals[1 * N.ncols + i] = fix[i].value[1];
		N.vals[2 * N.ncols + i] = fix[i].value[2];
	}
	getNtNNtN(N, NtNNtN[domain]);

	R.resize(3, K.nrows);
	int i = 0;
	for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++i) {
		for (auto di = dmap->begin(); di != dmap->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				switch (i % 2) {
				case 0:
					R.vals[0 * R.ncols + di->index] = 1;
					R.vals[1 * R.ncols + di->index] = 0;
					R.vals[2 * R.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 2].y;
					break;
				case 1:
					R.vals[0 * R.ncols + di->index] = 0;
					R.vals[1 * R.ncols + di->index] = 1;
					R.vals[2 * R.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 2].x;
					break;
				}
			}
		}
	}
}

// RedMat = Nt * ((N * Nt)^-1 * N)
//
//      per Fix-Point
// N = [  1  0  0
//        0  1  0
//        0  0  1
//       -y  x  0
//       -z  0  x
//        0 -z  y ]
//

template <typename T>
void RegularizationElasticity<T>::set3D(esint domain)
{
	std::vector<esint> fixPoints;
	this->getFixPoints(fixPoints, domain);

	const Matrix_CSR<T> &K = feti.K[domain];
	Matrix_Dense<T> &R = feti.R1[domain];
	Matrix_CSR<T> &RegMat = feti.RegMat[domain];
	const FETIDecomposition *decomposition = feti.decomposition;

	struct __fix__ { int col; double value[6]; }; // reorder DOFs to get sorted output
	std::vector<__fix__> fix(3 * fixPoints.size());
	for (size_t i = 0; i < fixPoints.size(); ++i) {
		auto dmap0 = decomposition->dmap->begin() + (3 * fixPoints[i] + 0);
		auto dmap1 = decomposition->dmap->begin() + (3 * fixPoints[i] + 1);
		auto dmap2 = decomposition->dmap->begin() + (3 * fixPoints[i] + 2);
		for (auto di = dmap0->begin(); di != dmap0->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				fix[3 * i + 0].col = di->index;
				fix[3 * i + 0].value[0] = 1;
				fix[3 * i + 0].value[1] = 0;
				fix[3 * i + 0].value[2] = 0;
				fix[3 * i + 0].value[3] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].y;
				fix[3 * i + 0].value[4] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].z;
				fix[3 * i + 0].value[5] = 0;
			}
		}
		for (auto di = dmap1->begin(); di != dmap1->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				fix[3 * i + 1].col = di->index;
				fix[3 * i + 1].value[0] = 0;
				fix[3 * i + 1].value[1] = 1;
				fix[3 * i + 1].value[2] = 0;
				fix[3 * i + 1].value[3] =  info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].x;
				fix[3 * i + 1].value[4] = 0;
				fix[3 * i + 1].value[5] = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].z;
			}
		}
		for (auto di = dmap2->begin(); di != dmap2->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				fix[3 * i + 2].col = di->index;
				fix[3 * i + 2].value[0] = 0;
				fix[3 * i + 2].value[1] = 0;
				fix[3 * i + 2].value[2] = 1;
				fix[3 * i + 2].value[3] = 0;
				fix[3 * i + 2].value[4] =  info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].x;
				fix[3 * i + 2].value[5] =  info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].y;
			}
		}
	}
	std::sort(fix.begin(), fix.end(), [] (const __fix__ &f1, const __fix__ &f2) { return f1.col < f2.col; });
	setRegMat<__fix__, T>(fix, RegMat, K.nrows);

	Matrix_Dense<T> N;
	N.resize(6, fix.size());
	for (size_t i = 0; i < fix.size(); ++i) {
		N.vals[0 * N.ncols + i] = fix[i].value[0];
		N.vals[1 * N.ncols + i] = fix[i].value[1];
		N.vals[2 * N.ncols + i] = fix[i].value[2];
		N.vals[3 * N.ncols + i] = fix[i].value[3];
		N.vals[4 * N.ncols + i] = fix[i].value[4];
		N.vals[5 * N.ncols + i] = fix[i].value[5];
	}
	getNtNNtN(N, NtNNtN[domain]);

	R.resize(6, K.nrows);
	int i = 0;
	for (auto dmap = decomposition->dmap->cbegin(); dmap != decomposition->dmap->cend(); ++dmap, ++i) {
		for (auto di = dmap->begin(); di != dmap->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				switch (i % 3) {
				case 0:
					R.vals[0 * R.ncols + di->index] = 1;
					R.vals[1 * R.ncols + di->index] = 0;
					R.vals[2 * R.ncols + di->index] = 0;
					R.vals[3 * R.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 3].y;
					R.vals[4 * R.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 3].z;
					R.vals[5 * R.ncols + di->index] = 0;
					break;
				case 1:
					R.vals[0 * R.ncols + di->index] = 0;
					R.vals[1 * R.ncols + di->index] = 1;
					R.vals[2 * R.ncols + di->index] = 0;
					R.vals[3 * R.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 3].x;
					R.vals[4 * R.ncols + di->index] = 0;
					R.vals[5 * R.ncols + di->index] = -info::mesh->nodes->coordinates->datatarray()[i / 3].z;
					break;
				case 2:
					R.vals[0 * R.ncols + di->index] = 0;
					R.vals[1 * R.ncols + di->index] = 0;
					R.vals[2 * R.ncols + di->index] = 1;
					R.vals[3 * R.ncols + di->index] = 0;
					R.vals[4 * R.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 3].x;
					R.vals[5 * R.ncols + di->index] =  info::mesh->nodes->coordinates->datatarray()[i / 3].y;
					break;
				}
			}
		}
	}
}

template <typename T>
void RegularizationElasticity<T>::setAnalytic()
{
	NtNNtN.resize(feti.K.size());
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		switch (info::mesh->dimension) {
		case 2: set2D(d); break;
		case 3: set3D(d); break;
		}
	}
}

template <typename T>
void RegularizationElasticity<T>::updateAnalytic()
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		const Matrix_CSR<T> &K = feti.K[d];
		Matrix_CSR<T> &RegMat = feti.RegMat[d];

		double max = 0;
		for (esint r = 0; r < K.nrows; ++r) {
			max = std::max(max, K.vals[K.rows[r] - Indexing::CSR]);
		}

		for (esint r = 0, i = 0; r < NtNNtN[d].nrows; ++r) {
			for (esint c = r; c < NtNNtN[d].ncols; ++c, ++i) {
				RegMat.vals[i] = max * NtNNtN[d].vals[r * NtNNtN[d].ncols + c];
			}
		}
	}
}

}
