
#include "fetisystem.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"

void espreso::initGluing(AX_FETISystem<double> &solver)
{
//	B1Dirichlet.initDomains(K.domains);
//	B1c.initDomains(DataDecomposition::DUPLICATION::SPLIT, K.domains);
//	B1Gluing.initDomains(K.domains);
//	B1Inequality.initDomains(K.domains);
//	B1duplication.initDomains(DataDecomposition::DUPLICATION::SPLIT, K.domains);
//	B1gap.initDomains(DataDecomposition::DUPLICATION::DUPLICATE, K.domains);
//
//	esint doffset = 0;
//	auto map = K.dmap->cbegin();
//	for (esint i = 0, prev = 0; i < BC.holder()->nnz; prev = BC.holder()->indices[i++]) {
//		map += BC.holder()->indices[i] - prev;
//		for (auto di = map->begin(); di != map->end(); ++di) {
//			if (K.ismy(di->domain)) { ++doffset; }
//		}
//	}
//
//	// dirichlet size (offset for the first gluing lambda)
//	esint dsize = Communication::exscan(doffset);
//
//	std::vector<std::vector<esint> > ROWS(K.domains), COLS(K.domains);
//	std::vector<std::vector<double> > VALS(K.domains), DUPS(K.domains), GAPS(K.domains);
//
//	map = K.dmap->cbegin();
//	for (esint i = 0, prev = 0; i < BC.holder()->nnz; prev = BC.holder()->indices[i++]) {
//		map += BC.holder()->indices[i] - prev;
//		for (auto di = map->begin(); di != map->end(); ++di) {
//			if (K.ismy(di->domain)) {
//				B1Map.push_back(doffset);
//				B1Map.push_back(0);
//				ROWS[di->domain - K.doffset].push_back(1 + doffset++);
//				COLS[di->domain - K.doffset].push_back(1 + di->index);
//				VALS[di->domain - K.doffset].push_back(1);
//			}
//		}
//	}
//
//	for (esint d = 0; d < K.domains; ++d) {
//		B1Dirichlet[d].resize(dsize, K[d].nrows, ROWS[d].size());
//		B1Dirichlet[d].fillPattern(ROWS[d].size(), ROWS[d].data(), COLS[d].data());
//		B1Dirichlet[d].fillValues(VALS[d].size(), VALS[d].data());
//		B1c[d].resize(VALS[d].size());
//		B1c[d].fill(0);
//
//		ROWS[d].clear();
//		COLS[d].clear();
//		VALS[d].clear();
//	}
//
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
//
//	esint goffset = 0;
//	std::vector<std::vector<esint> > sBuffer(info::mesh->neighbors.size()), rBuffer(info::mesh->neighbors.size());
//
//	auto fillbuffer = [&] (espreso::serializededata<esint, espreso::DI>::const_iterator &dmap) -> bool {
//		if (K.ismy(dmap->begin()->domain)) {
//			esint roffset = 0, prev = -1;
//			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
//				if (!K.ismy(di->domain)) {
//					while (K.distribution[K.neighbors[roffset] + 1] <= di->domain) {
//						++roffset;
//					}
//					if (prev != roffset) {
//						sBuffer[roffset].push_back(goffset);
//					}
//					prev = roffset;
//				}
//			}
//			return true;
//		}
//		return false;
//	};
//
//	size_t maxMultiplicity = 2;
//	if (K.nshared) {
//		map = K.dmap->cbegin();
//		for (esint n = 0, j = 0, prev = 0; n < K.nshared; prev = K.shared[n++]) {
//			map += K.shared[n] - prev;
//			while (j < BC.holder()->nnz && BC.holder()->indices[j] < K.shared[n]) {
//				++j;
//			}
//			if (j == BC.holder()->nnz || BC.holder()->indices[j] != K.shared[n]) {
//				maxMultiplicity = std::max(map->size(), maxMultiplicity);
//				if (fillbuffer(map)) {
//					if (solver.configuration.redundant_lagrange) {
//						goffset += map->size() * (map->size() - 1) / 2;
//					} else {
//						goffset += map->size() - 1;
//					}
//				}
//			}
//		}
//	}
//
//	esint gsize = Communication::exscan(goffset) + dsize + mortarRows;
//	goffset += dsize + mortarRows;
//
//	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
//		for (size_t i = 0; i < sBuffer[n].size(); ++i) {
//			sBuffer[n][i] += goffset;
//		}
//	}
//
//	if (!Communication::receiveLowerUnknownSize(sBuffer, rBuffer, info::mesh->neighbors)) {
//		eslog::internalFailure("exchange gluing offsets.\n");
//	}
//
//	MatrixDense multiplicity(maxMultiplicity - 1, maxMultiplicity);
//	for (esint r = 0; r < multiplicity.nrows; ++r) {
//		multiplicity[r][r] = 1;
//		multiplicity[r][r + 1] = -1;
//	}
//	multiplicity.orthonormalizeRows();
//
//	if (K.nshared) {
//		std::vector<esint> boffset(info::mesh->neighbors.size());
//		map = K.dmap->cbegin();
//		for (esint n = 0, j = 0, prev = 0; n < K.nshared; prev = K.shared[n++]) {
//			map += K.shared[n] - prev;
//			while (j < BC.holder()->nnz && BC.holder()->indices[j] < K.shared[n]) {
//				++j;
//			}
//			if (j == BC.holder()->nnz || BC.holder()->indices[j] != K.shared[n]) {
//				esint lambda;
//				if (K.ismy(map->begin()->domain)) {
//					lambda = goffset;
//					if (solver.configuration.redundant_lagrange) {
//						goffset += map->size() * (map->size() - 1) / 2;
//					} else {
//						goffset += map->size() - 1;
//					}
//				} else {
//					esint roffset = 0;
//					while (K.distribution[K.neighbors[roffset] + 1] <= map->begin()->domain) {
//						++roffset;
//					}
//					lambda = rBuffer[roffset][boffset[roffset]++];
//				}
//				if (solver.configuration.redundant_lagrange) {
//					for (auto di1 = map->begin(); di1 != map->end(); ++di1) {
//						for (auto di2 = di1 + 1; di2 != map->end(); ++di2) {
//							if (K.ismy(di1->domain) || K.ismy(di2->domain)) {
//								B1Map.push_back(lambda);
//								B1Map.push_back(0);
//
//								if (K.ismy(di1->domain)) {
//									ROWS[di1->domain - K.doffset].push_back(lambda + 1);
//									COLS[di1->domain - K.doffset].push_back(di1->index + 1);
//									VALS[di1->domain - K.doffset].push_back(1);
//									DUPS[di1->domain - K.doffset].push_back(1. / map->size());
//								} else {
//									esint roffset = 0;
//									while (K.distribution[K.neighbors[roffset] + 1] <= di1->domain) {
//										++roffset;
//									}
//									++B1Map.back();
//									B1Map.push_back(K.neighbors[roffset]);
//								}
//								if (K.ismy(di2->domain)) {
//									ROWS[di2->domain - K.doffset].push_back(lambda + 1);
//									COLS[di2->domain - K.doffset].push_back(di2->index + 1);
//									VALS[di2->domain - K.doffset].push_back(-1);
//									DUPS[di2->domain - K.doffset].push_back(1. / map->size());
//								} else {
//									esint roffset = 0;
//									while (K.distribution[K.neighbors[roffset] + 1] <= di2->domain) {
//										++roffset;
//									}
//									++B1Map.back();
//									B1Map.push_back(K.neighbors[roffset]);
//								}
//							}
//							++lambda;
//						}
//					}
//				} else {
//					// orthonormal multiplicity fill the lower triangle in the following way
//					// myfirst = the first dof that this process has
//					// x x 0 0 0
//					// x x x 0 0
//					// x x x x 0
//					// x x x x x
//					size_t myfirst = 0;
//					for (size_t i = 0; i < map->size(); ++i) {
//						if (K.ismy(map->at(i).domain)) {
//							myfirst = i;
//							break;
//						}
//					}
//					for (size_t link = 0; link + 1 < map->size(); ++link) {
//						if (myfirst <= link + 1) {
//							B1Map.push_back(lambda);
//							size_t mapCounter = B1Map.size();
//							B1Map.push_back(0);
//							esint roffset = 0;
//							for (size_t i = 0; i <= link + 1; ++i) {
//								if (K.ismy(map->at(i).domain)) {
//									ROWS[map->at(i).domain - K.doffset].push_back(lambda + 1);
//									COLS[map->at(i).domain - K.doffset].push_back(map->at(i).index + 1);
//									VALS[map->at(i).domain - K.doffset].push_back(multiplicity[link][i]);
//									DUPS[map->at(i).domain - K.doffset].push_back(1. / map->size());
//								} else {
//									while (K.distribution[K.neighbors[roffset] + 1] <= map->at(i).domain) {
//										++roffset;
//									}
//									if (B1Map[mapCounter] == 0 || B1Map.back() != K.neighbors[roffset]) {
//										++B1Map[mapCounter];
//										B1Map.push_back(K.neighbors[roffset]);
//									}
//								}
//							}
//						}
//						++lambda;
//					}
//				}
//			}
//		}
//	}
//
//	for (esint d = 0; d < K.domains; ++d) {
//		B1Gluing[d].resize(gsize, K[d].nrows, ROWS[d].size());
//		B1Gluing[d].fillPattern(ROWS[d].size(), ROWS[d].data(), COLS[d].data());
//		B1Gluing[d].fillValues(VALS[d].size(), VALS[d].data());
//		B1duplication[d].resize(DUPS[d].size());
//		B1duplication[d].fillValues(DUPS[d].data());
//
//		ROWS[d].clear();
//		COLS[d].clear();
//		VALS[d].clear();
//	}
//
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

void espreso::initKernels(AX_FETISystem<double> &solver, AX_HeatTransfer &assembler)
{
	solver.N1.domains.resize(solver.A.domains.size());
	solver.N2.domains.resize(solver.A.domains.size());
	solver.RegMat.domains.resize(solver.A.domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < solver.A.domains.size(); ++d) {
		if (assembler.hasKernel(d)) {
			solver.N1.domains[d].resize(solver.A.domains[d].nrows, 1);
			solver.N1.domains[d].type = Matrix_Type::REAL_UNSYMMETRIC;

			solver.RegMat.domains[d].resize(solver.A.domains[d].nrows, solver.A.domains[d].ncols, 1);
			solver.RegMat.domains[d].type = solver.A.domains[d].type;

			solver.RegMat.domains[d].rows[0] = solver.RegMat.domains[d].cols[0] = _Matrix_CSR_Pattern::Indexing;
			std::fill(solver.RegMat.domains[d].rows + 1, solver.RegMat.domains[d].rows + solver.RegMat.domains[d].nrows + 1, _Matrix_CSR_Pattern::Indexing + 1);
		}
	}
}

void espreso::updateKernels(AX_FETISystem<double> &solver, AX_HeatTransfer &assembler)
{
	#pragma omp parallel for
	for (size_t d = 0; d < solver.A.domains.size(); ++d) {
		if (assembler.hasKernel(d)) {
			Vector_Dense<double> diag;
			diag.resize(solver.A.domains[d].nrows);
			math::getDiagonal(solver.A.domains[d], diag);
			solver.RegMat.domains[d].vals[0] = math::max(diag);

			if (solver.configuration.method != FETIConfiguration::METHOD::HYBRID_FETI) { // N1 orthogonal for whole cluster
				esint crows = 0;
				for (esint dd = 0; dd < info::mesh->domains->size; dd++) {
					if (info::mesh->domains->cluster[d] == info::mesh->domains->cluster[dd]) {
						crows += solver.A.domains[dd].nrows;
					}
				}
				math::fill(solver.N1.domains[d], 1.0 / std::sqrt(crows));
			} else {
				math::fill(solver.N1.domains[d], 1.0 / std::sqrt(solver.A.domains[d].nrows));
			}
		}
	}
}
