
#include "regularization.elasticity.h"

#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/fetidatastore.h"
#include "wrappers/metis/w.metis.h"

#include <algorithm>
#include <numeric>
#include <random>

namespace espreso {

template struct RegularizationElasticity<double>;

static void getFixPoints(std::vector<esint> &fixPoints, int domain)
{
	esint begin = info::mesh->domains->elements[domain];
	esint end = info::mesh->domains->elements[domain + 1];

	esint FIX_POINTS_SIZE = 8;

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

	std::vector<std::vector<esint> > pids(FIX_POINTS_SIZE), pdist(FIX_POINTS_SIZE, { 1 }), pdata(FIX_POINTS_SIZE);
	for (size_t i = 0; i < partition.size(); i++) {
		pids[partition[i]].push_back(i);
	}
	for (size_t i = 0; i < partition.size(); i++) {
		esint p = partition[i];
		for (esint j = dist[i]; j < dist[i + 1]; j++) {
			if (partition[data[j]] == p) {
				size_t index = std::lower_bound(pids[p].begin(), pids[p].end(), data[j]) - pids[p].begin();
				if (pdist[p].size() <= index) {
					pdata[p].push_back(index + 1);
				}
			}
		}
		pdist[p].push_back(pdata[p].size() + 1);
	}

	for (esint p = 0; p < FIX_POINTS_SIZE; p++) {
		if (pids[p].size()) {
			std::vector<float> vals(pdata[p].size(), 1), x(pids[p].size(), 1. / pids[p].size()), y(pids[p].size());
			Matrix_CSR<float> M;
			M.shape = Matrix_Shape::UPPER;
			Vector_Dense<float> in, out;
			in.size = out.size = M.nrows = M.ncols = pids[p].size();
			M.nnz = vals.size();
			M.rows = pdist[p].data();
			M.cols = pdata[p].data();
			M.vals = vals.data();
			in.vals = x.data();
			out.vals = y.data();
			SpBLAS<float, Matrix_CSR> spblas(M);

			float last_l = pids[p].size(), l = 1;
			while (fabs((l - last_l) / l) > 1e-6) {
				spblas.apply(out, 1, 0, in);
				last_l = l;
				l = math::norm(out);
				math::scale(1 / l, out);
				out.swap(in);
			}
			fixPoints.push_back(ids[pids[p][std::max_element(in.vals, in.vals + in.size) - in.vals]]);
		}
	}
	utils::sortAndRemoveDuplicates(fixPoints);
}


// in  = Nt
// out = RegMat
//
// out = N * ((Nt * N)^-1 * Nt)
// out = N * (   A   )^-1 * Nt)
// out = N * (       X        )
template <typename T>
static void NNtNNt(Matrix_CSR<T> &in, Matrix_CSR<T> &out)
{
	SpBLAS<T, Matrix_CSR> A, X, N, Nt;
	Matrix_CSR<T> mA;

	Nt.insert(in);
	N.insertTransposed(in);
	A.multiply(Nt, N);
	A.extractUpper(mA);
}

template <typename T>
void RegularizationElasticity<T>::set2D(esint domain)
{
	std::vector<esint> fixPoints;
	getFixPoints(fixPoints, domain);

	const Matrix_CSR<T> &K = feti.K.domains[domain];
	Matrix_Dense<T> &R = feti.regularization.R1.domains[domain];
	Matrix_CSR<T> &RegMat = feti.regularization.RegMat.domains[domain];
	const DOFsDecomposition *decomposition = feti.K.decomposition;

	std::vector<std::pair<esint, T> > fixDofs(2 * fixPoints.size());
	for (size_t i = 0; i < fixPoints.size(); ++i) {
		auto dmap0 = decomposition->dmap->begin() + (2 * fixPoints[i] + 0);
		auto dmap1 = decomposition->dmap->begin() + (2 * fixPoints[i] + 1);
		for (auto di = dmap0->begin(); di != dmap0->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				fixDofs[i].first = -info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].y;
				fixDofs[i].second = di->index;
				break;
			}
		}
		for (auto di = dmap1->begin(); di != dmap1->end(); ++di) {
			if (di->domain == domain + decomposition->dbegin) {
				fixDofs[fixPoints.size() + i].first = info::mesh->nodes->coordinates->datatarray()[fixPoints[i]].x;
				fixDofs[fixPoints.size() + i].second = di->index;
				break;
			}
		}
	}
	std::sort(fixDofs.begin(), fixDofs.begin() + fixPoints.size(), [] (const std::pair<T, esint> &p1, const std::pair<T, esint> &p2) {
		return p1.second < p2.second;
	});
	std::sort(fixDofs.begin() + fixPoints.size(), fixDofs.end(), [] (const std::pair<T, esint> &p1, const std::pair<T, esint> &p2) {
		return p1.second < p2.second;
	});

	R.resize(K.nrows, 3);

	Matrix_CSR<T> Nt;
	Nt.resize(3, K.ncols, 4 * fixPoints.size());

	Nt.rows[0] = Indexing::CSR;
	Nt.rows[1] = Nt.rows[0] + fixPoints.size();
	Nt.rows[2] = Nt.rows[1] + fixPoints.size();
	Nt.rows[3] = Nt.rows[2] + 2 * fixPoints.size();

	esint cindex = 0;
	for (size_t i = 0; i < fixDofs.size(); i++, cindex++) {
		Nt.cols[cindex] = fixDofs[i].second + Indexing::CSR;
		Nt.vals[cindex] = 1;
	}

	std::sort(fixDofs.begin(), fixDofs.end(), [] (const std::pair<T, esint> &p1, const std::pair<T, esint> &p2) {
		return p1.second < p2.second;
	});

	for (size_t i = 0; i < fixDofs.size(); ++i, ++cindex) {
		Nt.cols[cindex] = fixDofs[i].second + Indexing::CSR;
		Nt.vals[cindex] = fixDofs[i].first;
	}

	NNtNNt(Nt, RegMat);

//	Nt.transposeTo(&N);
//	Matrix_CSR<T> A;
//	A.multiply(Nt, N);
//	A.removeLower(MatrixType::REAL_SYMMETRIC_INDEFINITE);
//	MatrixDense _X(N.nrows, N.ncols), B = N;
//	A.solve(B, _X);
//	MatrixCSR X = _X;
//	X.transpose();
//	_RegMat->at(d)->multiply(N, X);
//	_RegMat->at(d)->removeLower(MatrixType::REAL_SYMMETRIC_INDEFINITE);
//
//	RegMat[d].shallowCopyStructure(_RegMat->at(d));

	// Gram-schmidt
//	N1[d].resize(K[d].nrows, 3);
}

template <typename T>
void RegularizationElasticity<T>::set3D(esint domain)
{

}


template <typename T>
void RegularizationElasticity<T>::setAnalytic()
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		switch (info::mesh->dimension) {
		case 2: set2D(d); break;
		case 3: set3D(d); break;
		}
	}
}

template <typename T>
void RegularizationElasticity<T>::updateAnalytic()
{

}

}
