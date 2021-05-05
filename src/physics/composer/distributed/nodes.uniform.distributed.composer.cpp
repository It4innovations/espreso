
#include "nodes.uniform.distributed.composer.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"
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

NodesUniformDistributedComposer::NodesUniformDistributedComposer(Kernel *kernel, DistributedAssemblerData *data, int DOFs)
: DistributedComposer(kernel, data), _DOFs(DOFs)
{

}

void NodesUniformDistributedComposer::init()
{
	eslog::startln("DISTRIBUTED COMPOSER: STARTED", "COMPOSER");
	_initDOFMap();
	eslog::checkpoint("COMPOSER: DOF MAP COMPUTED");
	eslog::param("source", "NODE");
	eslog::param("DOFs/NODE", _DOFs);
	eslog::ln();

	_buildPatterns();
	clearMatrices(Builder::Request::KCM | Builder::Request::RBCf, _data);
	eslog::checkpointln("COMPOSER: MATRICES PATTERNS CREATED");

	_buildDirichlet();
	eslog::endln("COMPOSER: DIRICHLET INDICES COMPUTED");
}

void NodesUniformDistributedComposer::_initDOFMap()
{
	_DOFMap = new serializededata<esint, esint>(_DOFs, info::mesh->nodes->distribution);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		for (size_t n = info::mesh->nodes->distribution[t]; n < info::mesh->nodes->distribution[t + 1]; ++n) {
			for (int dof = 0; dof < _DOFs; ++dof) {
				_DOFMap->datatarray()[_DOFs * n + dof] = _DOFs * info::mesh->nodes->uniqInfo.position[n] + dof;
			}
		}
	}
}

void NodesUniformDistributedComposer::_buildPatterns()
{
	eslog::startln("COMPOSER: BUILD PATTERS", "PATTERNS BUILDING");
	size_t threads = info::env::OMP_NUM_THREADS;
	bool omitLower = false;

	_nDistribution = info::mesh->nodes->gatherUniqueNodeDistribution();
	for (size_t n = 0; n < _nDistribution.size(); ++n) {
		_nDistribution[n] *= _DOFs;
	}
	esint foreignDOFs = _DOFs * info::mesh->nodes->uniqInfo.nhalo;

	eslog::checkpointln("COMPOSER: GATHER DISTRIBUTION");

	_tKOffsets.resize(threads);
	_tRHSOffsets.resize(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		esint tKsize = 0, tRHSsize = 0;
		for (auto e = info::mesh->elements->nodes->begin(t); e != info::mesh->elements->nodes->end(t); ++e) {
			tRHSsize += e->size() * _DOFs;
			tKsize += getMatrixSize(e->size() * _DOFs, omitLower);
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (kernel->boundaryWithSettings(r)) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for (auto e = info::mesh->boundaryRegions[r]->elements->begin(t); e != info::mesh->boundaryRegions[r]->elements->end(t); ++e) {
						tRHSsize += e->size() * _DOFs;
						tKsize += getMatrixSize(e->size() * _DOFs, omitLower);
					}
				} else {
					for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(t); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(t); ++n) {
						if (info::mesh->nodes->size - info::mesh->nodes->uniqInfo.size <= *n) {
							tRHSsize += _DOFs;
						}
					}
				}
			}
		}

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
	esint localRHSOffset = utils::sizesToOffsets(_tRHSOffsets);

	eslog::checkpointln("COMPOSER: LOCAL OFFSETS");

	std::vector<IJ> KPattern(localKOffset);
	std::vector<esint> RHSPattern(localRHSOffset);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		IJ *Koffset = KPattern.data() + _tKOffsets[t];
		esint *RHSoffset = RHSPattern.data() + _tRHSOffsets[t];

		auto insert = [&] (serializededata<esint, esint>::const_iterator &e) {
			esint *_RHS = RHSoffset;
			for (esint dof = 0; dof < _DOFs; ++dof) {
				for (auto n = e->begin(); n != e->end(); ++n, ++RHSoffset) {
					*RHSoffset = _DOFMap->datatarray()[*n * _DOFs + dof];
				}
			}
			insertKPattern(Koffset, _RHS, RHSoffset, omitLower);
		};

		for (auto e = info::mesh->elements->nodes->cbegin(t); e != info::mesh->elements->nodes->cend(t); ++e) {
			insert(e);
			Koffset += getMatrixSize(e->size() * _DOFs, omitLower);
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->dimension && kernel->boundaryWithSettings(r)) {
				for (auto e = info::mesh->boundaryRegions[r]->elements->cbegin(t); e != info::mesh->boundaryRegions[r]->elements->cend(t); ++e) {
					insert(e);
					Koffset += getMatrixSize(e->size() * _DOFs, omitLower);
				}
			}
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
			if (info::mesh->boundaryRegions[r]->dimension == 0 && kernel->boundaryWithSettings(r)) {
				for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().cbegin(t); n != info::mesh->boundaryRegions[r]->nodes->datatarray().cend(t); ++n) {
					if (info::mesh->nodes->size - info::mesh->nodes->uniqInfo.size <= *n) {
						for (esint dof = 0; dof < _DOFs; ++dof, ++RHSoffset) {
							*RHSoffset = _DOFMap->datatarray()[*n * _DOFs + dof];
						}
					}
				}
			}
		}

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
		eslog::internalFailure("exchange K pattern.\n");
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

	_KPermutation.resize(KPattern.size());
	_RHSPermutation.resize(RHSPattern.size());

	std::vector<esint> ROW, COL;
	std::vector<double> VAL;

	ROW.reserve(info::mesh->nodes->size * _DOFs + 1);
	ROW.push_back(1);
	COL.push_back(KPattern[pK.front()].column + 1);
	_KPermutation[pK.front()] = 0;
	_RHSPermutation[pRHS.front()] = 0;
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
	for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
		if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
			++nonzeros;
		}
		_RHSPermutation[pRHS[i]] = nonzeros;
	}
	ROW.push_back(COL.size() + 1);
	VAL.resize(COL.size());

	eslog::checkpointln("COMPOSER: PATTERNS BUILD");

	_data->K.type = kernel->solverDataProvider->general->getMatrixType();
	_data->K.resize(info::mesh->nodes->size * _DOFs, info::mesh->nodes->uniqInfo.totalSize * _DOFs, COL.size(), foreignDOFs, info::mesh->neighbors.size());
	_data->K.fillPattern(_data->K.nrows, ROW.data(), COL.data());
	_data->K.fillDistribution(_DOFMap->datatarray().data(), _nDistribution.data(), info::mesh->neighbors.data());
	_data->K.structureUpdated();

	_data->M.shallowCopyStructure(&_data->K);
	_data->C.shallowCopyStructure(&_data->K);
	_data->CM.shallowCopyStructure(&_data->K);

	// TODO: share distribution with K
	_data->f.initVectors(kernel->solutions.size());
	_data->f.resize(info::mesh->nodes->size * _DOFs, foreignDOFs, info::mesh->neighbors.size());
	_data->f.fillDistribution(_DOFMap->datatarray().data(), _nDistribution.data(), info::mesh->neighbors.data());
	_data->f.structureUpdated();
	_data->R.shallowCopyStructure(&_data->f);
	_data->x.shallowCopyStructure(&_data->f);

	eslog::endln("COMPOSER: PATTERN FILLED");
}

void NodesUniformDistributedComposer::_buildDirichlet()
{
	std::vector<std::pair<esint, esint> > dIndices;
	kernel->solverDataProvider->general->dirichletIndices(dIndices);

	for (auto i = dIndices.begin(); i != dIndices.end(); ++i) {
		_dirichletMap.push_back(i->first * _DOFs + i->second);
	}

	_dirichletPermutation.resize(_dirichletMap.size());
	std::iota(_dirichletPermutation.begin(), _dirichletPermutation.end(), 0);
	std::sort(_dirichletPermutation.begin(), _dirichletPermutation.end(), [&] (esint i, esint j) {
		return _dirichletMap[i] < _dirichletMap[j];
	});

	std::sort(_dirichletMap.begin(), _dirichletMap.end());

	std::vector<esint> dmap = _dirichletMap;
	utils::removeDuplicates(dmap);

	_data->BC.initVectors(1);
	_data->BC.resize(info::mesh->nodes->size * _DOFs, dmap.size());
	_data->BC.holder()->fillPattern(dmap.data());
}
