
#include "analysis/linearsystem/fetisystem.h"
#include "basis/utilities/utils.h"
#include "esinfo/meshinfo.h"
#include "math2/generalization/matrix_feti.h"
#include "math2/primitives/matrix_dense.h"
#include "mesh/store/domainstore.h"

#include <algorithm>

namespace espreso {

template <typename T>
void _composeEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename AX_FETI<T>::EqualityConstraints &equalityConstraints, bool redundantLagrange)
{
	equalityConstraints.B1.domains.resize(K.domains.size());
	equalityConstraints.B1c.domains.resize(K.domains.size());
	equalityConstraints.B1Duplication.domains.resize(K.domains.size());
//	B1gap.initDomains(DataDecomposition::DUPLICATION::DUPLICATE, K.domains);
	equalityConstraints.D2C.resize(K.domains.size());

	std::vector<esint> lsize, l2mpi;

	esint doffset = 0;
	auto map = K.decomposition->dmap->cbegin();
	for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
		map += dirichlet.cluster.indices[i] - prev;
		for (auto di = map->begin(); di != map->end(); ++di) {
			if (K.decomposition->ismy(di->domain)) { ++doffset; }
		}
	}

	// dirichlet size (offset for the first gluing lambda)
	esint dsize = Communication::exscan(doffset), coffset = 0;

	std::vector<std::vector<esint> > ROWS(K.domains.size()), COLS(K.domains.size());
	std::vector<std::vector<T> > VALS(K.domains.size()), DUPS(K.domains.size());

	map = K.decomposition->dmap->cbegin();
	for (esint i = 0, prev = 0; i < dirichlet.cluster.nnz; prev = dirichlet.cluster.indices[i++]) {
		map += dirichlet.cluster.indices[i] - prev;
		for (auto di = map->begin(); di != map->end(); ++di) {
			if (K.decomposition->ismy(di->domain)) {
				lsize.push_back(0);
				equalityConstraints.C2G.push_back(doffset++);
				equalityConstraints.D2C[di->domain - K.decomposition->dbegin].push_back(coffset++);
				ROWS[di->domain - K.decomposition->dbegin].push_back(1 + ROWS[di->domain - K.decomposition->dbegin].size());
				COLS[di->domain - K.decomposition->dbegin].push_back(1 + di->index);
				VALS[di->domain - K.decomposition->dbegin].push_back(1);
			}
		}
	}

//	for (size_t d = 0; d < K.domains.size(); ++d) {
//		equalityConstraints.B1Dirichlet.domains[d].resize(dsize, K.domains[d].nrows, ROWS[d].size());
//		equalityConstraints.B1c.domains[d].resize(VALS[d].size());
//		std::copy(ROWS[d].begin(), ROWS[d].end(), equalityConstraints.B1Dirichlet.domains[d].rows);
//		std::copy(COLS[d].begin(), COLS[d].end(), equalityConstraints.B1Dirichlet.domains[d].cols);
//		std::copy(VALS[d].begin(), VALS[d].end(), equalityConstraints.B1Dirichlet.domains[d].vals);
//		std::fill(equalityConstraints.B1c.domains[d].vals, equalityConstraints.B1c.domains[d].vals + equalityConstraints.B1c.domains[d].size, 0);
//
//		ROWS[d].clear();
//		COLS[d].clear();
//		VALS[d].clear();
//	}

//	// MORTAR matrices
//	esint mortarRows = 0;
//	std::vector<esint> mlambdas;
//	for (esint d = 0; d < mortars.domains; ++d) {
//		for (esint i = 0; i < mortars[d].nnz; ++i) {
//			ROWS[d].push_back(mortars[d].rows[i] + dsize);
//			mlambdas.push_back(mortars[d].rows[i] + dsize - 1);
//			if (i == 0 || mortars[d].rows[i - 1] != mortars[d].rows[i]) {
//				DUPS[d].push_back(0.5);
//			}
//		}
//		COLS[d].insert(COLS[d].end(), mortars[d].cols, mortars[d].cols + mortars[d].nnz);
//		VALS[d].insert(VALS[d].end(), mortars[d].vals, mortars[d].vals + mortars[d].nnz);
//	}
//
//	if (mortars.domains) {
//		mortarRows = mortars[0].nrows;
//		utils::sortAndRemoveDuplicates(mlambdas);
//		auto dmap = mortars.dmap->begin();
//		for (size_t i = 0; i < mlambdas.size(); ++i, ++dmap) {
//			B1Map.push_back(mlambdas[i]);
//			size_t size = B1Map.size();
//			B1Map.push_back(0);
//			esint roffset = 0, prev = -1;
//			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
//				if (!mortars.ismy(di->domain)) {
//					while (mortars.distribution[mortars.neighbors[roffset] + 1] <= di->domain) {
//						++roffset;
//					}
//					if (prev != roffset) {
//						++B1Map[size];
//						B1Map.push_back(mortars.neighbors[roffset]);
//					}
//					prev = roffset;
//				}
//			}
//		}
//	}
//
////	Communication::serialize([&] () {
////		std::cout << info::mpi::rank << "\n";
////		std::cout << ROWS[0];
////		std::cout << COLS[0];
////		std::cout << VALS[0];
////		std::cout << "MAP: " << B1Map;
////	});

	esint goffset = 0, lastNeigh = info::mesh->neighbors.size();
	std::vector<std::vector<esint> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());

	auto fillbuffer = [&] (espreso::serializededata<esint, espreso::DIndex>::const_iterator &dmap) -> bool {
		if (K.decomposition->ismy(dmap->begin()->domain)) {
			esint roffset = 0, prev = -1;
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (!K.decomposition->ismy(di->domain)) {
					while (roffset + 1 < lastNeigh && K.decomposition->neighDomain[roffset + 1] <= di->domain) {
						++roffset;
					}
					if (prev != roffset) {
						sBuffer[roffset].push_back(goffset);
					}
					prev = roffset;
				}
			}
			return true;
		}
		return false;
	};

	size_t maxMultiplicity = 2;
	if (K.decomposition->sharedDOFs.size()) {
		map = K.decomposition->dmap->cbegin();
		esint dindex = 0;
		for (size_t n = 0, prev = 0; n < K.decomposition->sharedDOFs.size(); prev = K.decomposition->sharedDOFs[n++]) {
			map += K.decomposition->sharedDOFs[n] - prev;
			while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < K.decomposition->sharedDOFs[n]) {
				++dindex;
			}
			if (dindex == dirichlet.cluster.nnz || dirichlet.cluster.indices[dindex] != K.decomposition->sharedDOFs[n]) {
				maxMultiplicity = std::max(map->size(), maxMultiplicity);
				if (fillbuffer(map)) {
					if (redundantLagrange) {
						goffset += map->size() * (map->size() - 1) / 2;
					} else {
						goffset += map->size() - 1;
					}
				}
			}
		}
	}

	esint gsize = Communication::exscan(goffset) + dsize;
	goffset += dsize;

	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		for (size_t i = 0; i < sBuffer[n].size(); ++i) {
			sBuffer[n][i] += goffset;
		}
	}

	if (!Communication::receiveLowerUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
		eslog::internalFailure("exchange gluing offsets.\n");
	}

	Matrix_Dense<double> multiplicity;
	multiplicity.resize(maxMultiplicity - 1, maxMultiplicity);
	math::set(multiplicity, 0.0);
	for (esint r = 0; r < multiplicity.nrows; ++r) {
		multiplicity.vals[r * multiplicity.ncols + r] = 1;
		multiplicity.vals[r * multiplicity.ncols + r + 1] = -1;
	}
	math::orthonormalize(multiplicity);

	if (K.decomposition->sharedDOFs.size()) {
		std::vector<esint> boffset(info::mesh->neighbors.size());
		map = K.decomposition->dmap->cbegin();
		esint dindex = 0;
		for (size_t n = 0, prev = 0; n < K.decomposition->sharedDOFs.size(); prev = K.decomposition->sharedDOFs[n++]) {
			map += K.decomposition->sharedDOFs[n] - prev;
			while (dindex < dirichlet.cluster.nnz && dirichlet.cluster.indices[dindex] < K.decomposition->sharedDOFs[n]) {
				++dindex;
			}
			if (dindex == dirichlet.cluster.nnz || dirichlet.cluster.indices[dindex] != K.decomposition->sharedDOFs[n]) {
				esint lambda;
				if (K.decomposition->ismy(map->begin()->domain)) {
					lambda = goffset;
					if (redundantLagrange) {
						goffset += map->size() * (map->size() - 1) / 2;
					} else {
						goffset += map->size() - 1;
					}
				} else {
					esint roffset = 0;
					while (roffset + 1 < lastNeigh && K.decomposition->neighDomain[roffset + 1] <= map->begin()->domain) {
						++roffset;
					}
					lambda = rBuffer[roffset][boffset[roffset]++];
				}
				if (redundantLagrange) {
					for (auto di1 = map->begin(); di1 != map->end(); ++di1) {
						for (auto di2 = di1 + 1; di2 != map->end(); ++di2) {
							if (K.decomposition->ismy(di1->domain) || K.decomposition->ismy(di2->domain)) {
								equalityConstraints.C2G.push_back(lambda);
								lsize.push_back(0);

								if (K.decomposition->ismy(di1->domain)) {
									equalityConstraints.D2C[di1->domain - K.decomposition->dbegin].push_back(coffset);
									ROWS[di1->domain - K.decomposition->dbegin].push_back(1 + ROWS[di1->domain - K.decomposition->dbegin].size());
									COLS[di1->domain - K.decomposition->dbegin].push_back(di1->index + 1);
									VALS[di1->domain - K.decomposition->dbegin].push_back(1);
									DUPS[di1->domain - K.decomposition->dbegin].push_back(1. / map->size());
								} else {
									esint roffset = 0;
									while (roffset + 1 < lastNeigh && K.decomposition->neighDomain[roffset + 1] <= di1->domain) {
										++roffset;
									}
									++lsize.back();
									l2mpi.push_back(roffset);
								}
								if (K.decomposition->ismy(di2->domain)) {
									equalityConstraints.D2C[di2->domain - K.decomposition->dbegin].push_back(coffset);
									ROWS[di2->domain - K.decomposition->dbegin].push_back(1 + ROWS[di2->domain - K.decomposition->dbegin].size());
									COLS[di2->domain - K.decomposition->dbegin].push_back(di2->index + 1);
									VALS[di2->domain - K.decomposition->dbegin].push_back(-1);
									DUPS[di2->domain - K.decomposition->dbegin].push_back(1. / map->size());
								} else {
									esint roffset = 0;
									while (roffset + 1 < lastNeigh && K.decomposition->neighDomain[roffset + 1] <= di2->domain) {
										++roffset;
									}
									++lsize.back();
									l2mpi.push_back(roffset);
								}
								++coffset;
							}
							++lambda;
						}
					}
				} else {
					// orthonormal multiplicity fill the lower triangle in the following way
					// myfirst = the first dof that this process has
					// x x 0 0 0
					// x x x 0 0
					// x x x x 0
					// x x x x x
					size_t myfirst = 0;
					for (size_t i = 0; i < map->size(); ++i) {
						if (K.decomposition->ismy(map->at(i).domain)) {
							myfirst = i;
							break;
						}
					}
					for (size_t link = 0; link + 1 < map->size(); ++link) {
						if (myfirst <= link + 1) {
							equalityConstraints.C2G.push_back(lambda);
							lsize.push_back(0);
							for (size_t i = 0; i <= link + 1; ++i) {
								if (K.decomposition->ismy(map->at(i).domain)) {
									equalityConstraints.D2C[map->at(i).domain - K.decomposition->dbegin].push_back(coffset);
									ROWS[map->at(i).domain - K.decomposition->dbegin].push_back(1 + ROWS[map->at(i).domain - K.decomposition->dbegin].size());
									COLS[map->at(i).domain - K.decomposition->dbegin].push_back(map->at(i).index + 1);
									VALS[map->at(i).domain - K.decomposition->dbegin].push_back(multiplicity.vals[link * multiplicity.ncols + i]);
									DUPS[map->at(i).domain - K.decomposition->dbegin].push_back(1. / map->size());
								} else {
									esint roffset = 0;
									while (roffset + 1 < lastNeigh && K.decomposition->neighDomain[roffset + 1] <= map->at(i).domain) {
										++roffset;
									}
									if (lsize.back() == 0 || l2mpi.back() != roffset) {
										++lsize.back();
										l2mpi.push_back(roffset);
									}
								}
							}
						}
						++lambda;
					}
				}
			}
		}
	}

	for (size_t d = 0; d < K.domains.size(); ++d) {
		equalityConstraints.B1.domains[d].resize(gsize, K.domains[d].nrows, ROWS[d].size());
		equalityConstraints.B1Duplication.domains[d].resize(DUPS[d].size());

		std::copy(ROWS[d].begin(), ROWS[d].end(), equalityConstraints.B1.domains[d].rows);
		std::copy(COLS[d].begin(), COLS[d].end(), equalityConstraints.B1.domains[d].cols);
		std::copy(VALS[d].begin(), VALS[d].end(), equalityConstraints.B1.domains[d].vals);
		std::copy(DUPS[d].begin(), DUPS[d].end(), equalityConstraints.B1Duplication.domains[d].vals);

		ROWS[d].clear();
		COLS[d].clear();
		VALS[d].clear();
	}

	lsize.push_back(0);
	utils::sizesToOffsets(lsize);
	equalityConstraints.L2MPI = new serializededata<esint, int>(lsize, l2mpi);

//	esint ioffset = 0;
//	map = K.dmap->cbegin();
//	for (esint i = 0, prev = 0; i < gapDirection.holder()->nnz; prev = gapDirection.holder()->indices[i++]) {
//		map += gapDirection.holder()->indices[i] - prev;
//		if (K.ismy(map->begin()->domain)) { ++ioffset; }
//	}
//
//	esint isize = Communication::exscan(ioffset) + gsize;
//	ioffset += gsize;
//
//	map = K.dmap->cbegin();
//	for (esint i = 0, prev = 0; i < gapDirection.holder()->nnz; prev = gapDirection.holder()->indices[i++]) {
//		map += gapDirection.holder()->indices[i] - prev;
//		if (K.ismy(map->begin()->domain)) {
//			B1Map.push_back(ioffset);
//			B1Map.push_back(0);
//			ROWS[map->begin()->domain - K.doffset].push_back(1 + ioffset++);
//			COLS[map->begin()->domain - K.doffset].push_back(1 + map->begin()->index);
//			VALS[map->begin()->domain - K.doffset].push_back(gapDirection.holder()->vals[i]);
//			GAPS[map->begin()->domain - K.doffset].push_back(gap.holder()->vals[i]);
//		}
//	}
//
//	for (esint d = 0; d < K.domains; ++d) {
//		B1Inequality[d].resize(isize, K[d].nrows, ROWS[d].size());
//		B1Inequality[d].fillPattern(ROWS[d].size(), ROWS[d].data(), COLS[d].data());
//		B1Inequality[d].fillValues(VALS[d].size(), VALS[d].data());
//		B1gap[d].resize(GAPS[d].size());
//		B1gap[d].fillValues(GAPS[d].data());
//
//		ROWS[d].clear();
//		COLS[d].clear();
//		VALS[d].clear();
//		GAPS[d].clear();
//	}
}

template <typename T>
void _evaluateEqualityConstraints(const Matrix_FETI<Matrix_CSR, T> &K, const Vector_Distributed<Vector_Sparse, T> &dirichlet, typename AX_FETI<T>::EqualityConstraints &equalityConstraints, bool redundantLagrange)
{

}

void composeHeatTransferKernel(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
{
	R.resize(K.nrows, 1);
	R.type = Matrix_Type::REAL_NONSYMMETRIC;
	R.shape = Matrix_Shape::FULL;

	RegMat.resize(K.nrows, K.ncols, 1);
	RegMat.type = K.type;
	RegMat.shape = K.shape;

	RegMat.rows[0] = RegMat.cols[0] = _Matrix_CSR_Pattern::Indexing;
	std::fill(RegMat.rows + 1, RegMat.rows + RegMat.nrows + 1, RegMat.rows[0] + 1);
}

void evaluateHeatTransferKernel(const Matrix_CSR<double> &K, Matrix_Dense<double> &R, Matrix_CSR<double> &RegMat)
{
	RegMat.vals[0] = math::getDiagonalMax(K);
	math::set(R, 1.0 / std::sqrt(K.nrows));

//	if (solver.configuration.method != FETIConfiguration::METHOD::HYBRID_FETI) { // N1 orthogonal for whole cluster
//		esint crows = 0;
//		for (esint dd = 0; dd < info::mesh->domains->size; dd++) {
//			if (info::mesh->domains->cluster[d] == info::mesh->domains->cluster[dd]) {
//				crows += solver.A.domains[dd].nrows;
//			}
//		}
//		math::fill(solver.N1.domains[d], 1.0 / std::sqrt(crows));
//	} else {
//		math::fill(solver.N1.domains[d], 1.0 / std::sqrt(solver.A.domains[d].nrows));
//	}
}

template <> void composeEqualityConstraints(const Matrix_FETI<Matrix_CSR, double> &K, const Vector_Distributed<Vector_Sparse, double> &dirichlet, AX_FETI<double>::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _composeEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }
template <> void composeEqualityConstraints(const Matrix_FETI<Matrix_CSR, std::complex<double> > &K, const Vector_Distributed<Vector_Sparse, std::complex<double> > &dirichlet, AX_FETI<std::complex<double> >::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _composeEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }

template <> void evaluateEqualityConstraints(const Matrix_FETI<Matrix_CSR, double> &K, const Vector_Distributed<Vector_Sparse, double> &dirichlet, AX_FETI<double>::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _evaluateEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }
template <> void evaluateEqualityConstraints(const Matrix_FETI<Matrix_CSR, std::complex<double> > &K, const Vector_Distributed<Vector_Sparse, std::complex<double> > &dirichlet, AX_FETI<std::complex<double> >::EqualityConstraints &equalityConstraints, bool redundantLagrange) { _evaluateEqualityConstraints(K, dirichlet, equalityConstraints, redundantLagrange); }

}
