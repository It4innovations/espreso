
#include "regularization.elasticity.h"

#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "math/wrappers/math.blas.h"
#include "math/wrappers/math.spblas.h"
#include "math/wrappers/math.lapack.h"
#include "wrappers/metis/w.metis.h"

#include "math/math.h"

#include <algorithm>
#include <numeric>
#include <random>

namespace espreso {

template <typename T>
void RegularizationElasticity<T>::getFixPoints(std::vector<esint> &fixPoints, int domain)
{
	esint begin = info::mesh->domains->elements[domain];
	esint end = info::mesh->domains->elements[domain + 1];

	size_t FIX_POINTS_SIZE = 8;

	auto neighs = [] (std::vector<esint> &neighs, Element::CODE code, int node, const esint* nodes) {
		switch (code) {
		case Element::CODE::HEXA8:
		case Element::CODE::HEXA20:
			if (node < 4) {
				neighs.push_back((nodes[(node + 1) % 4]));
				neighs.push_back((nodes[(node + 3) % 4]));
				neighs.push_back((nodes[node + 4]));
			} else {
				neighs.push_back((nodes[(node + 1) % 4 + 4]));
				neighs.push_back((nodes[(node + 3) % 4 + 4]));
				neighs.push_back((nodes[node - 4]));
			}
			return 3;
		case Element::CODE::TETRA4:
		case Element::CODE::TETRA10:
			neighs.push_back(nodes[(node + 1) % 4]);
			neighs.push_back(nodes[(node + 2) % 4]);
			neighs.push_back(nodes[(node + 3) % 4]);
			return 3;
		case Element::CODE::PRISMA6:
		case Element::CODE::PRISMA15:
			if (node < 3) {
				neighs.push_back(nodes[(node + 1) % 3]);
				neighs.push_back(nodes[(node + 2) % 3]);
				neighs.push_back(nodes[node + 3]);
			} else {
				neighs.push_back(nodes[(node + 1) % 3 + 3]);
				neighs.push_back(nodes[(node + 2) % 3 + 3]);
				neighs.push_back(nodes[node - 3]);
			}
			return 3;

		case Element::CODE::PYRAMID5:
		case Element::CODE::PYRAMID13:
			if (node == 4) {
				neighs.insert(neighs.end(), nodes, nodes + 4);
				return 4;
			} else {
				neighs.push_back(nodes[(node + 1) % 4]);
				neighs.push_back(nodes[(node + 3) % 4]);
				neighs.push_back(nodes[4]);
				return 3;
			}

		case Element::CODE::TRIANGLE3:
		case Element::CODE::TRIANGLE6:
			neighs.push_back(nodes[(node + 1) % 3]);
			neighs.push_back(nodes[(node + 2) % 3]);
			return 2;

		case Element::CODE::SQUARE4:
		case Element::CODE::SQUARE8:
			neighs.push_back(nodes[(node + 1) % 4]);
			neighs.push_back(nodes[(node + 3) % 4]);
			return 2;

		case Element::CODE::LINE2:
		case Element::CODE::LINE3:
			neighs.push_back(nodes[(node + 1) % 2]);
			return 1;
		case Element::CODE::POINT1:
		default:
			return 0;
		}
		return 0;
	};

	std::vector<esint> originnodes, neighsnodes;
	originnodes.reserve((end - begin) * 20);
	auto element = info::mesh->elements->nodes->begin() + begin;
	const auto &epointer = info::mesh->elements->epointers->datatarray();
	for (esint e = 0; e < end - begin; ++e, ++element) {
		for (int n = 0; n < epointer[begin + e]->coarseNodes; ++n) {
			originnodes.insert(
					originnodes.end(),
					neighs(neighsnodes, epointer[begin + e]->code, n, element->data()),
					element->at(n));
		}
	}

	std::vector<esint> permutation(originnodes.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		return originnodes[i] < originnodes[j];
	});

	std::vector<esint> ids, dist, data;
	dist.push_back(0);
	ids.push_back(originnodes[permutation[0]]);
	data.push_back(neighsnodes[permutation[0]]);
	for (size_t i = 1; i < permutation.size(); i++) {
		if (originnodes[permutation[i]] != originnodes[permutation[i - 1]]) {
			utils::sortAndRemoveDuplicates(data, dist.back());
			dist.push_back(data.size());
			ids.push_back(originnodes[permutation[i]]);
		}
		data.push_back(neighsnodes[permutation[i]]);
	}
	utils::sortAndRemoveDuplicates(data, dist.back());
	dist.push_back(data.size());

	if (ids.size() <= FIX_POINTS_SIZE * FIX_POINTS_SIZE / 2) {
		std::random_device rd;
		std::mt19937 g(rd());

		std::shuffle(ids.begin(), ids.end(), g);
		if (FIX_POINTS_SIZE < ids.size()) {
			ids.resize(FIX_POINTS_SIZE);
		}
		fixPoints = ids;
		utils::sortAndRemoveDuplicates(fixPoints);
		return;
	}

	for (size_t i = 0; i < data.size(); i++) {
		data[i] = std::lower_bound(ids.begin(), ids.end(), data[i]) - ids.begin();
	}

	std::vector<esint> partition(ids.size());
	METISConfiguration options;
	METIS::call(options, ids.size(), dist.data(), data.data(), 0, NULL, NULL, FIX_POINTS_SIZE, partition.data());

	std::vector<std::vector<esint> > pids(FIX_POINTS_SIZE), pdist(FIX_POINTS_SIZE, { Indexing::CSR }), pdata(FIX_POINTS_SIZE);
	for (size_t i = 0; i < partition.size(); i++) {
		pids[partition[i]].push_back(i);
	}
	for (size_t i = 0; i < partition.size(); i++) {
		esint p = partition[i];
		for (esint j = dist[i]; j < dist[i + 1]; j++) {
			if (partition[data[j]] == p) {
				size_t index = std::lower_bound(pids[p].begin(), pids[p].end(), data[j]) - pids[p].begin();
				if (pdist[p].size() <= index) {
					pdata[p].push_back(index + Indexing::CSR);
				}
			}
		}
		pdist[p].push_back(pdata[p].size() + Indexing::CSR);
	}

	for (size_t p = 0; p < FIX_POINTS_SIZE; p++) {
		if (pids[p].size()) {
			std::vector<float> vals(pdata[p].size(), 1), x(pids[p].size(), 1. / pids[p].size()), y(pids[p].size());
			Matrix_CSR<float> M;
			M.shape = Matrix_Shape::UPPER;
			M.type = Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE;
			Vector_Dense<float> in, out;
			in.size = out.size = M.nrows = M.ncols = pids[p].size();
			M.nnz = vals.size();
			M.rows = pdist[p].data();
			M.cols = pdata[p].data();
			M.vals = vals.data();
			in.vals = x.data();
			out.vals = y.data();
			SpBLAS<Matrix_CSR, float> spblas(M);

			float last_l = pids[p].size(), l = 1;
			while (fabs((l - last_l) / l) > 1e-6) {
				spblas.apply(out, 1, 0, in);
				last_l = l;
				l = math::blas::norm(out.size, out.vals, 1);
				math::blas::scale(out.size, 1 / l, out.vals, 1);
				out.swap(in);
			}
			fixPoints.push_back(ids[pids[p][std::max_element(in.vals, in.vals + in.size) - in.vals]]);
		}
	}
	utils::sortAndRemoveDuplicates(fixPoints);
}

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
void RegularizationElasticity<T>::set2D(FETI<T> &feti, esint domain)
{
	std::vector<esint> fixPoints;
	getFixPoints(fixPoints, domain);

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
void RegularizationElasticity<T>::set3D(FETI<T> &feti, esint domain)
{
	std::vector<esint> fixPoints;
	getFixPoints(fixPoints, domain);

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
void RegularizationElasticity<T>::set(FETI<T> &feti)
{
	NtNNtN.resize(feti.K.size());
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.size(); ++d) {
		switch (info::mesh->dimension) {
		case 2: set2D(feti, d); break;
		case 3: set3D(feti, d); break;
		}
	}
}

template <typename T>
void RegularizationElasticity<T>::update(FETI<T> &feti)
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

template <typename T>
std::vector<Matrix_Dense<T> > RegularizationElasticity<T>::NtNNtN;

template struct RegularizationElasticity<double>;

}
