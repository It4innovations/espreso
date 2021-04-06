
#include "faces.edges.uniform.distributed.composer.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/communication.h"
#include "basis/containers/serializededata.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "physics/kernels/kernel.h"
#include "physics/kernels/solverdataprovider/provider.h"
#include "physics/system/distributedsystem.h"
#include "physics/system/builder/builder.h"
#include "math/matrix.indices.h"
#include "math/matrix.csr.distributed.h"
#include "math/vector.sparse.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

FacesEdgesUniformDistributedComposer::FacesEdgesUniformDistributedComposer(Kernel *kernel, DistributedAssemblerData *data, int fDOFs, int eDOFs)
: DistributedComposer(kernel, data), _fDOFs(fDOFs), _eDOFs(eDOFs)
{

}

void FacesEdgesUniformDistributedComposer::init()
{
	_initDOFMap();
	eslog::checkpoint("COMPOSER: DATA MAPPED TO DOFs");
	eslog::param("source", "FACE, EDGE");
	eslog::param("DOFs/FACE", _fDOFs);
	eslog::param("DOFs/EDGE", _eDOFs);
	eslog::ln();

	_buildPatterns();
	clearMatrices(Builder::Request::KCM | Builder::Request::RBCf, _data);
	eslog::checkpointln("COMPOSER: MATRICES PATTERNS CREATED");

	_buildDirichlet();
	eslog::checkpointln("COMPOSER: DIRICHLET INDICES COMPUTED");
}

void FacesEdgesUniformDistributedComposer::_initDOFMap()
{
	int threads = info::env::OMP_NUM_THREADS;

	std::vector<esint> faceOffset(threads), edgeOffset(threads);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		esint dsize = 0;
		auto eID = info::mesh->elements->IDs->datatarray().cbegin(t);
		auto faces = info::mesh->elements->faceNeighbors->cbegin(t);
		for (size_t n = info::mesh->elements->distribution[t]; n < info::mesh->elements->distribution[t + 1]; ++n, ++faces, ++eID) {
			for (auto face = faces->begin(); face != faces->end(); ++face) {
				if (*face == -1 || *eID < *face) {
					dsize += _fDOFs;
				}
			}
		}
		faceOffset[t] = dsize;

		dsize = 0;
		eID = info::mesh->elements->IDs->datatarray().cbegin(t);
		auto edges = info::mesh->elements->edgeNeighbors->cbegin(t);
		for (size_t n = info::mesh->elements->distribution[t]; n < info::mesh->elements->distribution[t + 1]; ++n, ++edges, ++eID) {
			for (auto edge = edges->begin(); edge != edges->end(); edge += *edge + 1) {
				if (*edge == 0 || *eID < *(edge + 1)) {
					dsize += _eDOFs;
				}
			}
		}
		edgeOffset[t] = dsize;
	}

	esint foffset = utils::sizesToOffsets(faceOffset);
	esint goffset = foffset + utils::sizesToOffsets(edgeOffset);
	for (int t = 0; t < threads; ++t) {
		edgeOffset[t] += foffset;
	}
	_nDistribution = Communication::getDistribution(goffset);
	goffset = _nDistribution[info::mpi::rank];

	std::vector<esint> edistribution = Communication::getDistribution(info::mesh->elements->size);

	std::vector<std::vector<std::vector<esint> > > sIDBuffer(threads, std::vector<std::vector<esint> >(info::mesh->neighbors.size()));
	std::vector<std::vector<std::vector<esint> > > sDistBuffer(threads, std::vector<std::vector<esint> >(info::mesh->neighbors.size(), std::vector<esint>({0})));
	std::vector<std::vector<std::vector<esint> > > sDataBuffer(threads, std::vector<std::vector<esint> >(info::mesh->neighbors.size()));
	std::vector<std::vector<esint> > rIDBuffer(info::mesh->neighbors.size());
	std::vector<std::vector<esint> > rDistBuffer(info::mesh->neighbors.size());
	std::vector<std::vector<esint> > rDataBuffer(info::mesh->neighbors.size());

	std::vector<std::vector<esint> > DOFDist(threads), DOFData(threads);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		esint tfoffset = faceOffset[t];
		esint teoffset = edgeOffset[t];
		std::vector<esint> tdist, tdata;
		if (t == 0) {
			tdist.push_back(0);
		}

		auto eID = info::mesh->elements->IDs->datatarray().cbegin(t);
		auto faces = info::mesh->elements->faceNeighbors->cbegin(t);
		auto edges = info::mesh->elements->edgeNeighbors->cbegin(t);
		auto nodes = info::mesh->elements->procNodes->cbegin(t);
		for (size_t n = info::mesh->elements->distribution[t]; n < info::mesh->elements->distribution[t + 1]; ++n, ++faces, ++edges, ++nodes, ++eID) {
			for (auto face = faces->begin(); face != faces->end(); ++face) {
				if (*face != -1 && (*face < edistribution[info::mpi::rank] || edistribution[info::mpi::rank + 1] <= *face)) {
					esint noffset = 0;
					while (edistribution[info::mesh->neighbors[noffset] + 1] <= *face) {
						++noffset;
					}
					if (sIDBuffer[t][noffset].size() == 0 || sIDBuffer[t][noffset].back() != *eID) {
						sIDBuffer[t][noffset].push_back(*eID);
						sDistBuffer[t][noffset].push_back(sDistBuffer[t][noffset].back() + 2); // face offset
						sDistBuffer[t][noffset].push_back(sDistBuffer[t][noffset].back()); // edge offset
					} else {
						*(sDistBuffer[t][noffset].end() - 2) += 2;
						*(sDistBuffer[t][noffset].end() - 1) += 2;
					}
					sDataBuffer[t][noffset].push_back(*face);
					sDataBuffer[t][noffset].push_back(tfoffset);
				}
				if (*face == -1 || *eID < *face) {
					for (int dof = 0; dof < _fDOFs; ++dof) {
						tdata.push_back(goffset + tfoffset++);
					}
				} else {
					for (int dof = 0; dof < _fDOFs; ++dof) {
						tdata.push_back(-1); // push dummy
					}
				}
			}
			auto edgenodes = info::mesh->elements->epointers->datatarray()[n]->edges->cbegin();
			for (auto edge = edges->begin(); edge != edges->end(); edge += *edge + 1, ++edgenodes) {
				for (esint i = 1, prev = -1; i <= *edge; ++i) {
					if (edge[i] < edistribution[info::mpi::rank] || edistribution[info::mpi::rank + 1] <= edge[i]) {
						esint noffset = 0;
						while (edistribution[info::mesh->neighbors[noffset] + 1] <= edge[i]) {
							++noffset;
						}
						if (prev != noffset) {
							if (sIDBuffer[t][noffset].size() == 0 || sIDBuffer[t][noffset].back() != *eID) {
								sIDBuffer[t][noffset].push_back(*eID);
								sDistBuffer[t][noffset].push_back(sDistBuffer[t][noffset].back());
								sDistBuffer[t][noffset].push_back(sDistBuffer[t][noffset].back() + 3); // n1, n2, dof
							} else {
								sDistBuffer[t][noffset].back() += 3;
							}

							esint n1 = info::mesh->nodes->IDs->datatarray()[nodes->at(edgenodes->at(0))];
							esint n2 = info::mesh->nodes->IDs->datatarray()[nodes->at(edgenodes->at(1))];
							if (n1 < n2) {
								sDataBuffer[t][noffset].push_back(n1);
								sDataBuffer[t][noffset].push_back(n2);
							} else {
								sDataBuffer[t][noffset].push_back(n2);
								sDataBuffer[t][noffset].push_back(n1);
							}
							sDataBuffer[t][noffset].push_back(teoffset);
						}
						prev = noffset;
					}
				}
				if (*edge == 0 || *eID < *(edge + 1)) {
					for (int dof = 0; dof < _eDOFs; ++dof) {
						tdata.push_back(goffset + teoffset++);
					}
				} else {
					for (int dof = 0; dof < _eDOFs; ++dof) {
						tdata.push_back(-1); // push dummy
					}
				}
			}

			tdist.push_back(tdata.size());
		}
		DOFDist[t].swap(tdist);
		DOFData[t].swap(tdata);
	}

	utils::threadDistributionToFullDistribution(DOFDist);
	_DOFMap = new serializededata<esint, esint>(DOFDist, DOFData); // create map with empty spaces -> filled by the code belows

	for (int t = 1; t < threads; t++) {
		for (size_t n = 0; n < sIDBuffer[t].size(); n++) {
			sIDBuffer[0][n].insert(sIDBuffer[0][n].end(), sIDBuffer[t][n].begin(), sIDBuffer[t][n].end());
		}
		for (size_t n = 0; n < sDistBuffer[t].size(); n++) {
			sDistBuffer[0][n].insert(sDistBuffer[0][n].end(), sDistBuffer[t][n].begin(), sDistBuffer[t][n].end());
		}
		for (size_t n = 0; n < sDataBuffer[t].size(); n++) {
			sDataBuffer[0][n].insert(sDataBuffer[0][n].end(), sDataBuffer[t][n].begin(), sDataBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sIDBuffer[0], rIDBuffer, info::mesh->neighbors)) {
		eslog::error("ESPRESO internal error: exchange elements ID in faces-edges-distributed-composer.\n");
	}
	if (!Communication::exchangeUnknownSize(sDistBuffer[0], rDistBuffer, info::mesh->neighbors)) {
		eslog::error("ESPRESO internal error: exchange elements Dist in faces-edges-distributed-composer.\n");
	}
	if (!Communication::exchangeUnknownSize(sDataBuffer[0], rDataBuffer, info::mesh->neighbors)) {
		eslog::error("ESPRESO internal error: exchange elements Data in faces-edges-distributed-composer.\n");
	}

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		auto eID = info::mesh->elements->IDs->datatarray().cbegin(t);
		auto faces = info::mesh->elements->faceNeighbors->cbegin(t);
		auto edges = info::mesh->elements->edgeNeighbors->cbegin(t);
		auto nodes = info::mesh->elements->procNodes->cbegin(t);
		auto map = _DOFMap->datatarray().begin(t);
		for (size_t n = info::mesh->elements->distribution[t]; n < info::mesh->elements->distribution[t + 1]; ++n, ++faces, ++edges, ++nodes, ++eID) {
			for (auto face = faces->begin(); face != faces->end(); ++face, map += _fDOFs) {
				if (*face != -1 && *face < *eID) {
					if (*face < edistribution[info::mpi::rank] || edistribution[info::mpi::rank + 1] <= *face) {
						esint noffset = 0;
						while (edistribution[info::mesh->neighbors[noffset] + 1] <= *face) {
							++noffset;
						}
						size_t position = std::lower_bound(rIDBuffer[noffset].begin(), rIDBuffer[noffset].end(), *face) - rIDBuffer[noffset].begin();
						for (esint i = rDistBuffer[noffset][2 * position]; i < rDistBuffer[noffset][2 * position + 1]; i += 2) {
							if (rDataBuffer[noffset][i] == *eID) {
								for (int dof = 0; dof < _fDOFs; ++dof) {
									map[dof] = rDataBuffer[noffset][i + 1] + dof;
								}
								break;
							}
						}
					} else {
						auto nface = info::mesh->elements->faceNeighbors->cbegin() + *face;
						auto nmap = (_DOFMap->cbegin() + *face)->begin();
						for (auto nf = nface->begin(); nf != nface->end(); ++nf, nmap += _fDOFs) {
							if (*nf == *eID) {
								for (int dof = 0; dof < _fDOFs; ++dof) {
									map[dof] = nmap[dof];
								}
								break;
							}
						}
					}
				}
			}
			auto edgenodes = info::mesh->elements->epointers->datatarray()[n]->edges->cbegin();
			for (auto edge = edges->begin(); edge != edges->end(); edge += *edge + 1, ++edgenodes, map += _eDOFs) {
				if (edge[0] != 0 && edge[1] < *eID) {
					esint n1 = info::mesh->nodes->IDs->datatarray()[nodes->at(edgenodes->at(0))];
					esint n2 = info::mesh->nodes->IDs->datatarray()[nodes->at(edgenodes->at(1))];
					if (n2 < n1) { std::swap(n1, n2); }
					if (edge[1] < edistribution[info::mpi::rank] || edistribution[info::mpi::rank + 1] <= edge[1]) {
						esint noffset = 0;
						while (edistribution[info::mesh->neighbors[noffset] + 1] <= edge[1]) {
							++noffset;
						}
						size_t position = std::lower_bound(rIDBuffer[noffset].begin(), rIDBuffer[noffset].end(), edge[1]) - rIDBuffer[noffset].begin();
						for (esint i = rDistBuffer[noffset][2 * position + 1]; i < rDistBuffer[noffset][2 * position + 2]; i += 3) {
							if (rDataBuffer[noffset][i] == n1 && rDataBuffer[noffset][i + 1] == n2) {
								for (int dof = 0; dof < _fDOFs; ++dof) {
									map[dof] = rDataBuffer[noffset][i + 2] + dof;
								}
								break;
							}
						}
					} else {
						esint nfaces = (info::mesh->elements->faceNeighbors->cbegin() + edge[1])->size();
						auto nedge = info::mesh->elements->edgeNeighbors->cbegin() + edge[1];
						auto nnodes = info::mesh->elements->procNodes->cbegin() + edge[1];
						auto nedgenodes = (info::mesh->elements->epointers->datatarray()[edge[1]])->edges->cbegin();
						auto nmap = (_DOFMap->cbegin() + edge[1])->begin() + _fDOFs * nfaces;
						for (auto ne = nedge->begin(); ne != nedge->end(); ne += *ne + 1, ++nedgenodes, nmap += _eDOFs) {
							esint nn1 = info::mesh->nodes->IDs->datatarray()[nnodes->at(nedgenodes->at(0))];
							esint nn2 = info::mesh->nodes->IDs->datatarray()[nnodes->at(nedgenodes->at(1))];
							if (nn2 < nn1) { std::swap(nn1, nn2); }
							if (nn1 == n1 && nn2 == n2) {
								for (int dof = 0; dof < _eDOFs; ++dof) {
									map[dof] = nmap[dof];
								}
								break;
							}
						}
					}
				}
			}
		}
	}
}

void FacesEdgesUniformDistributedComposer::_buildPatterns()
{
	eslog::startln("COMPOSER: BUILD PATTERS", "PATTERNS BUILDING");
	int threads = info::env::OMP_NUM_THREADS;
	bool omitLower = false;

	eslog::checkpointln("COMPOSER: GATHER DISTRIBUTION");

	_tKOffsets.resize(threads);
	_tRHSOffsets.resize(threads);

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		esint tKsize = 0, tRHSsize = 0;

		auto map = _DOFMap->cbegin(t);
		for (size_t n = info::mesh->elements->distribution[t]; n < info::mesh->elements->distribution[t + 1]; ++n, ++map) {
			tRHSsize += map->size();
			tKsize += getMatrixSize(map->size(), omitLower);
		}

//		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
//			if (info::mesh->boundaryRegions[r]->dimension) {
//				for (auto e = info::mesh->boundaryRegions[r]->procNodes->begin(t); e != info::mesh->boundaryRegions[r]->procNodes->end(t); ++e) {
//					tRHSsize += e->size() * _DOFs;
//					tKsize += getMatrixSize(e->size() * _DOFs, omitLower);
//				}
//			} else {
//				for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(t); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(t); ++n) {
//					if (info::mesh->nodes->size - info::mesh->nodes->uniqueSize <= *n) {
//						tRHSsize += _DOFs;
//					}
//				}
//			}
//		}

//		for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
//			for (auto n = info::mesh->elementsRegions[r]->nodes->datatarray().begin(t); n != info::mesh->elementsRegions[r]->nodes->datatarray().end(t); ++n) {
//				if (info::mesh->nodes->size - info::mesh->nodes->uniqueSize <= *n) {
//					tRHSsize += _DOFs;
//				}
//			}
//		}

		_tKOffsets[t] = tKsize;
		_tRHSOffsets[t] = tRHSsize;
	}

	esint localKOffset = utils::sizesToOffsets(_tKOffsets);
	utils::sizesToOffsets(_tRHSOffsets);

	eslog::checkpointln("COMPOSER: LOCAL OFFSETS");

	std::vector<IJ> KPattern(localKOffset);
	std::vector<esint> RHSPattern(_DOFMap->datatarray().begin(), _DOFMap->datatarray().end());

	#pragma omp parallel for
	for (int t = 0; t < threads; t++) {
		IJ *Koffset = KPattern.data() + _tKOffsets[t];
		auto map = _DOFMap->cbegin(t);
		for (size_t n = info::mesh->elements->distribution[t]; n < info::mesh->elements->distribution[t + 1]; ++n, ++map) {
			insertKPattern(Koffset, map->begin(), map->end(), omitLower);
			Koffset += getMatrixSize(map->size(), omitLower);
		}

//		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
//			if (info::mesh->boundaryRegions[r]->dimension) {
//				for (auto e = info::mesh->boundaryRegions[r]->procNodes->cbegin(t); e != info::mesh->boundaryRegions[r]->procNodes->cend(t); ++e) {
//					insert(e);
//					Koffset += getMatrixSize(e->size() * _DOFs, omitLower);
//				}
//			}
//		}
//
//		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
//			if (info::mesh->boundaryRegions[r]->dimension == 0) {
//				for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().cbegin(t); n != info::mesh->boundaryRegions[r]->nodes->datatarray().cend(t); ++n) {
//					if (info::mesh->nodes->size - info::mesh->nodes->uniqueSize <= *n) {
//						for (esint dof = 0; dof < _DOFs; ++dof, ++RHSoffset) {
//							*RHSoffset = _DOFMap->datatarray()[*n * _DOFs + dof];
//						}
//					}
//				}
//			}
//		}

//		for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
//			for (auto n = info::mesh->elementsRegions[r]->nodes->datatarray().begin(t); n != info::mesh->elementsRegions[r]->nodes->datatarray().end(t); ++n) {
//				if (info::mesh->nodes->size - info::mesh->nodes->uniqueSize <= *n) {
//					for (esint dof = 0; dof < _DOFs; ++dof, ++RHSoffset) {
//						*RHSoffset = _DOFMap->datatarray()[*n * _DOFs + dof];
//					}
//				}
//			}
//		}
	}

	eslog::checkpointln("COMPOSER: COMPUTE LOCAL OFFSETS");

	std::vector<esint> pK(KPattern.size());
	std::iota(pK.begin(), pK.end(), 0);
	std::sort(pK.begin(), pK.end(), [&] (esint i, esint j) {
		return KPattern[i] < KPattern[j];
	});

	std::vector<esint> pRHS(RHSPattern.size());
	std::iota(pRHS.begin(), pRHS.end(), 0);
	std::sort(pRHS.begin(), pRHS.end(), [&] (esint i, esint j) {
		return RHSPattern[i] < RHSPattern[j];
	});

	eslog::checkpointln("COMPOSER: COMPUTE LOCAL PERMUTATION");

	std::vector<std::vector<IJ> > sKBuffer(info::mesh->neighbors.size()), rKBuffer(info::mesh->neighbors.size());

	auto iK = pK.begin();
	for (size_t n = 0; n < info::mesh->neighbors.size() && info::mesh->neighbors[n] < info::mpi::rank; ++n) {
		while (KPattern[*iK].row < _nDistribution[info::mesh->neighbors[n] + 1]) {
			if (iK == pK.begin() || KPattern[*iK] != KPattern[*(iK - 1)]) {
				sKBuffer[n].push_back(KPattern[*iK]);
			}
			++iK;
		}
	}

	eslog::checkpointln("COMPOSER: FILL SEND/RECV BUFFERS");

	if (!Communication::receiveUpperUnknownSize(sKBuffer, rKBuffer, info::mesh->neighbors)) {
		eslog::error("ESPRESO internal error: exchange K pattern.\n");
	}

	for (size_t i = 0; i < rKBuffer.size(); i++) {
		KPattern.insert(KPattern.end(), rKBuffer[i].begin(), rKBuffer[i].end());
	}

	eslog::checkpointln("COMPOSER: EXCHANGE BUFFERS");

	size_t localK = pK.size();
	pK.resize(KPattern.size());
	std::iota(pK.begin() + localK, pK.end(), localK);
	std::sort(pK.begin(), pK.end(), [&] (esint i, esint j) {
		return KPattern[i] < KPattern[j];
	});

	eslog::checkpointln("COMPOSER: COMPUTE NEIGHBORS PERMUTATION");

	_RHSPermutation.resize(RHSPattern.size());
	_RHSPermutation[pRHS.front()] = 0;
	size_t dofs = 0;
	for (size_t i = 1; i < pRHS.size(); i++) {
		if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
			++dofs;
		}
		_RHSPermutation[pRHS[i]] = dofs;
	}
	++dofs;

	_KPermutation.resize(KPattern.size());
	_KPermutation[pK.front()] = 0;

	std::vector<esint> ROW, COL;
	std::vector<double> VAL;
	ROW.reserve(dofs + 1);
	ROW.push_back(1);
	COL.push_back(KPattern[pK.front()].column + 1);
	for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
		if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
			++nonzeros;
			COL.push_back(KPattern[pK[i]].column + 1);
			if (KPattern[pK[i - 1]].row != KPattern[pK[i]].row) {
				ROW.push_back(nonzeros + 1);
			}
		}
		_KPermutation[pK[i]] = nonzeros;
	}

	ROW.push_back(COL.size() + 1);
	VAL.resize(COL.size());

	eslog::checkpointln("COMPOSER: PATTERNS BUILD");

	esint foreignDOFs = dofs - (_nDistribution[info::mpi::rank + 1] - _nDistribution[info::mpi::rank]);
	std::vector<esint> halo;
	halo.reserve(foreignDOFs);
	for (esint i = 0; i < foreignDOFs; i++) {
		halo.push_back(RHSPattern[pRHS[i]]);
	}

	_data->K.type = kernel->solverDataProvider->general->getMatrixType();
	_data->K.resize(dofs, _nDistribution.back(), COL.size(), foreignDOFs, info::mesh->neighbors.size());
	_data->K.fillPattern(_data->K.nrows, ROW.data(), COL.data());
	_data->K.fillDistribution(halo.data(), _nDistribution.data(), info::mesh->neighbors.data());
	_data->K.structureUpdated();

	_data->M.shallowCopyStructure(&_data->K);
	_data->C.shallowCopyStructure(&_data->K);

	// TODO: share distribution with K
	_data->f.initVectors(kernel->solutions.size());
	_data->f.resize(dofs, foreignDOFs, info::mesh->neighbors.size());
	_data->f.fillDistribution(halo.data(), _nDistribution.data(), info::mesh->neighbors.data());
	_data->f.structureUpdated();
	_data->R.shallowCopyStructure(&_data->f);
	_data->x.shallowCopyStructure(&_data->f);

	eslog::endln("COMPOSER: PATTERN FILLED");
}

void FacesEdgesUniformDistributedComposer::_buildDirichlet()
{
	std::vector<std::pair<esint, esint> > dIndices;
	kernel->solverDataProvider->general->dirichletIndices(dIndices);

//	for (auto i = dIndices.begin(); i != dIndices.end(); ++i) {
//		_dirichletMap.push_back(i->first * _DOFs + i->second);
//	}

	_dirichletPermutation.resize(_dirichletMap.size());
	std::iota(_dirichletPermutation.begin(), _dirichletPermutation.end(), 0);
	std::sort(_dirichletPermutation.begin(), _dirichletPermutation.end(), [&] (esint i, esint j) {
		return _dirichletMap[i] < _dirichletMap[j];
	});

	std::sort(_dirichletMap.begin(), _dirichletMap.end());

	std::vector<esint> dmap = _dirichletMap;
	utils::removeDuplicates(dmap);

	_data->BC.initVectors(1);
	_data->BC.resize(_data->f[0].size, dmap.size());
	_data->BC.holder()->fillPattern(dmap.data());
}
