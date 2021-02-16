
#include "nodes.uniform.feti.composer.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
#include "basis/utilities/utils.h"
#include "wrappers/mpi/communication.h"
#include "basis/containers/serializededata.h"
#include "config/ecf/input/decomposition.h"
#include "config/ecf/physics/physicssolver/loadstep.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/fetidatastore.h"
#include "mesh/store/contactstore.h"
#include "physics/kernels/kernel.h"
#include "physics/kernels/solverdataprovider/provider.h"
#include "physics/system/fetisystem.h"
#include "math/matrix.indices.h"
#include "math/domainindices.h"
#include "math/vector.sparse.h"
#include "output/visualization/debug.h"
#include "wrappers/metis/w.metis.h"

#include <algorithm>
#include <numeric>

using namespace espreso;

NodesUniformFETIComposer::NodesUniformFETIComposer(const FETIConfiguration &configuration, Kernel *kernel, FETIAssemblerData *data, int DOFs)
: FETIComposer(configuration, kernel, data), _DOFs(DOFs)
{

}

void NodesUniformFETIComposer::init()
{
	eslog::startln("FETI COMPOSER: STARTED", "COMPOSER");
	_initDOFMap();
	eslog::checkpoint("COMPOSER: NODES MAPPED TO DOFs");
	eslog::param("source", "NODE");
	eslog::param("DOFs/NODE", _DOFs);
	eslog::ln();

	if (_configuration.regularization == FETIConfiguration::REGULARIZATION::ANALYTIC) {
		computeFixPoints();
		eslog::checkpointln("COMPOSER: FIX POINTS COMPUTED");
	}

	if (_configuration.B0_type == FETIConfiguration::B0_TYPE::CORNERS) {
		computeCornerNodes();
		eslog::checkpointln("COMPOSER: CORNERS COMPUTED");
	}

	_buildPatterns();
	clearMatrices(Builder::Request::KCM | Builder::Request::RBCf, _data);
	eslog::checkpointln("COMPOSER: STRUCTURAL MATRICES CREATED");

	_buildDirichlet();
	eslog::checkpointln("COMPOSER: DIRICHLET INDICES COMPUTED");

	_buildMortars();
	eslog::checkpointln("COMPOSER: MORTART MATRICES INITIALIZED");

	_buildInequality();
	eslog::endln("COMPOSER: INEQUALITY CONSTRAINS COMPUTED");
}

void NodesUniformFETIComposer::_initDOFMap()
{
	// TODO: improve performance using info::mesh->nodes->domains
	size_t threads = info::env::OMP_NUM_THREADS;

	_domainDOFsSize.resize(info::mesh->elements->ndomains);

	// nID, domain
	std::vector<std::vector<std::pair<esint, esint> > > ntodomains(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<esint, esint> > tdata;

		for (esint d = info::mesh->elements->domainDistribution[t]; d != info::mesh->elements->domainDistribution[t + 1]; ++d) {
			std::vector<esint> dnodes(
					(info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d])->begin(),
					(info::mesh->elements->procNodes->begin() + info::mesh->elements->elementsDistribution[d + 1])->begin());

			utils::sortAndRemoveDuplicates(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				tdata.push_back(std::pair<esint, esint>(dnodes[i], d));
			}
			_domainDOFsSize[d] = dnodes.size() * _DOFs;
		}

		ntodomains[t].swap(tdata);
	}

	for (size_t t = 1; t < threads; t++) {
		ntodomains[0].insert(ntodomains[0].end(), ntodomains[t].begin(), ntodomains[t].end());
	}

	std::sort(ntodomains[0].begin(), ntodomains[0].end());

	std::vector<std::vector<esint> > DOFs(threads);
	std::vector<size_t> ndistribution = tarray<size_t>::distribute(threads, ntodomains[0].size());

	for (size_t t = 1; t < threads; t++) {
		while (
				ndistribution[t] < ntodomains[0].size() &&
				ntodomains[0][ndistribution[t]].first == ntodomains[0][ndistribution[t] - 1].first) {

			++ndistribution[t];
		}
	}

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> tDOFs(info::mesh->elements->ndomains);

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			tDOFs[ntodomains[0][n].second] += _DOFs;
		}

		DOFs[t].swap(tDOFs);
	}

	utils::sizesToOffsets(DOFs);

	std::vector<std::vector<std::vector<esint> > > sBuffer(threads);
	std::vector<std::vector<esint> > rBuffer(info::mesh->neighborsWithMe.size());

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto nranks = info::mesh->nodes->ranks->begin() + ntodomains[0][ndistribution[t]].first;

		std::vector<std::vector<esint> > tBuffer(info::mesh->neighborsWithMe.size());

		size_t n = ndistribution[t];
		while (n < ndistribution[t + 1]) {
			size_t begin = n++;
			while (n < ntodomains[0].size() && ntodomains[0][n].first == ntodomains[0][n - 1].first) {
				++n;
			}

			esint noffset = 0;
			for (auto r = nranks->begin(); r != nranks->end(); ++r) {
				while (info::mesh->neighborsWithMe[noffset] < *r) {
					++noffset;
				}

				tBuffer[noffset].push_back(n - begin);
				for (size_t i = begin; i < n; i++) {
					tBuffer[noffset].push_back(info::mesh->elements->firstDomain + ntodomains[0][i].second);
					tBuffer[noffset].push_back(DOFs[t][ntodomains[0][i].second]);
				}
			}
			++nranks;

			for (size_t i = begin; i < n; i++) {
				DOFs[t][ntodomains[0][i].second] += _DOFs;
			}
		}
		sBuffer[t].swap(tBuffer);
	}

	for (size_t t = 1; t < threads; t++) {
		for (size_t n = 0; n < sBuffer[t].size(); n++) {
			sBuffer[0][n].insert(sBuffer[0][n].end(), sBuffer[t][n].begin(), sBuffer[t][n].end());
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer[0], rBuffer, info::mesh->neighborsWithMe)) {
		eslog::internalFailure("exchange uniform DOFs.\n");
	}

	std::vector<esint> DOFDistribution(1);
	std::vector<DI> DOFData;

	// TODO: make it parallel
	// parallelization is possible if node order will be kept as: boundary first!
	// now we prefer generality
	auto nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> roffset(rBuffer.size());
	std::vector<DI> current; current.reserve(50);
	for (esint n = 0; n < info::mesh->nodes->size; ++n, ++nranks) {
		esint noffset = 0;
		for (auto r = nranks->begin(); r != nranks->end(); ++r) {
			while (info::mesh->neighborsWithMe[noffset] < *r) {
				++noffset;
			}

			esint domains = rBuffer[noffset][roffset[noffset]++];
			for (esint d = 0; d < domains; ++d) {
				current.push_back({rBuffer[noffset][roffset[noffset] + 2 * d], rBuffer[noffset][roffset[noffset] + 2 * d + 1]});
			}
			roffset[noffset] += 2 * domains;
		}

		for (int dof = 0; dof < _DOFs; ++dof) {
			for (size_t i = 0; i < current.size(); ++i) {
				DOFData.push_back({ current[i].domain, current[i].index + dof });
			}
			DOFDistribution.push_back(DOFData.size());
		}
		current.clear();
	}

	std::vector<size_t> distribution = info::mesh->nodes->distribution, datadistribution(threads + 1);
	for (size_t t = 1; t < threads; t++) {
		distribution[t] = _DOFs * distribution[t] + 1;
		datadistribution[t] = DOFDistribution[distribution[t]];
	}
	datadistribution[threads] = DOFDistribution[_DOFs * distribution[threads]];
	distribution[threads] = _DOFs * distribution[threads] + 1;

	_DOFMap = new serializededata<esint, DI>(
			tarray<esint>(distribution, DOFDistribution),
			tarray<DI>(datadistribution, DOFData));
}

void NodesUniformFETIComposer::_buildPatterns()
{
	_KPermutation.resize(info::mesh->elements->ndomains);
	_RHSPermutation.resize(info::mesh->elements->ndomains);

	std::vector<esint> distribution = info::mesh->elements->gatherDomainsProcDistribution();

	_data->K.initDomains(info::mesh->elements->ndomains);
	_data->K.fillDecomposition(
				info::mpi::rank, info::mpi::size, info::mesh->neighbors.size(),
				distribution.data(), info::mesh->neighbors.data(), _DOFMap);
	_data->f.initVectors(kernel->solutions.size());
	_data->f.initDomains(DataDecomposition::DUPLICATION::SPLIT, info::mesh->elements->ndomains);
	_data->f.fillDecomposition(
			info::mpi::rank, info::mpi::size, info::mesh->neighbors.size(),
			distribution.data(), info::mesh->neighbors.data(), _DOFMap);

	#pragma omp parallel for
	for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
		for (esint d = info::mesh->elements->domainDistribution[t]; d != info::mesh->elements->domainDistribution[t + 1]; ++d) {
			if (_BEMDomain[d]) {
				_buildKBEMPattern(d);
			} else {
				_buildKFEMPattern(d);
			}
			_data->f.resizeDomain(d, _domainDOFsSize[d]);
		}
	}

	_data->M.shallowCopyStructure(&_data->K);
	_data->C.shallowCopyStructure(&_data->K);
	_data->R.shallowCopyStructure(&_data->f);
	_data->x.shallowCopyStructure(&_data->f);
	_data->x.setDuplications(DataDecomposition::DUPLICATION::DUPLICATE);
}

void NodesUniformFETIComposer::_buildKFEMPattern(esint domain)
{
	_data->K[domain].type = kernel->solverDataProvider->feti->getMatrixType(domain);

	auto ebegin = info::mesh->elements->procNodes->cbegin() + info::mesh->elements->elementsDistribution[domain];
	auto eend = info::mesh->elements->procNodes->cbegin() + info::mesh->elements->elementsDistribution[domain + 1];

	esint Ksize = 0, RHSsize = 0;
	for (auto e = ebegin; e != eend; ++e) {
		RHSsize += e->size() * _DOFs;
		Ksize += getMatrixSize(e->size() * _DOFs, _data->K[domain].type != MatrixType::REAL_UNSYMMETRIC);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension && kernel->boundaryWithSettings(r)) {
			if (info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
				esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
				esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
				auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;

				for (esint i = begin; i < end; ++i, ++enodes) {
					RHSsize += enodes->size() * _DOFs;
					Ksize += getMatrixSize(enodes->size() * _DOFs, _data->K[domain].type != MatrixType::REAL_UNSYMMETRIC);
				}
			}
		}
	}
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (!info::mesh->boundaryRegions[r]->dimension && kernel->boundaryWithSettings(r)) {
			esint prev = 0;
			auto dmap = _DOFMap->begin();
			for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); prev = *n++) {
				dmap += (*n - prev) * _DOFs;
				if (dmap->begin()->domain == domain + info::mesh->elements->firstDomain) {
					RHSsize += _DOFs;
				}
			}
		}
	}

	std::vector<esint> permK(Ksize), permRHS(RHSsize);
	std::vector<IJ> KPattern(Ksize);
	std::vector<esint> RHSPattern(RHSsize), ROW, COL;

	IJ *Koffset = KPattern.data();
	esint *RHSoffset = RHSPattern.data();

	auto insert = [&] (serializededata<esint, esint>::const_iterator &enodes) {
		esint *_RHS = RHSoffset;
		for (int dof = 0; dof < _DOFs; ++dof) {
			for (auto n = enodes->begin(); n != enodes->end(); ++n, ++RHSoffset) {
				auto DOFs = (_DOFMap->begin() + (*n * _DOFs + dof))->begin();
				while (DOFs->domain != domain + info::mesh->elements->firstDomain) {
					++DOFs;
				}
				*RHSoffset = DOFs->index;
			}
		}
		insertKPattern(Koffset, _RHS, RHSoffset, _data->K[domain].type != MatrixType::REAL_UNSYMMETRIC);
	};

	for (auto e = ebegin; e != eend; ++e) {
		insert(e);
		Koffset += getMatrixSize(e->size() * _DOFs, _data->K[domain].type != MatrixType::REAL_UNSYMMETRIC);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension && kernel->boundaryWithSettings(r)) {
			if (info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
				esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
				esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
				auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;
				for (esint i = begin; i < end; ++i, ++enodes) {
					insert(enodes);
					Koffset += getMatrixSize(enodes->size() * _DOFs, _data->K[domain].type != MatrixType::REAL_UNSYMMETRIC);
				}
			}
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (!info::mesh->boundaryRegions[r]->dimension && kernel->boundaryWithSettings(r)) {
			esint prev = 0;
			auto dmap = _DOFMap->begin();
			for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); prev = *n++) {
				dmap += (*n - prev) * _DOFs;
				for (int dof = 0; dof < _DOFs; ++dof) {
					if ((dmap + dof)->begin()->domain == domain + info::mesh->elements->firstDomain) {
						*RHSoffset++ = (dmap + dof)->begin()->index;
					}
				}
			}
		}
	}

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

	ROW.reserve(_domainDOFsSize[domain] + 1);
	COL.reserve(KPattern.size());
	ROW.push_back(1);
	COL.push_back(KPattern[pK.front()].column + 1);
	permK[pK.front()] = 0;
	for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
		if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
			++nonzeros;
			COL.push_back(KPattern[pK[i]].column + 1);
			if (KPattern[pK[i - 1]].row != KPattern[pK[i]].row) {
				ROW.push_back(nonzeros + 1);
			}
		}
		permK[pK[i]] = nonzeros;
	}
	permRHS[pRHS.front()] = 0;
	for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
		if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
			++nonzeros;
		}
		permRHS[pRHS[i]] = nonzeros;
	}

	_KPermutation[domain].swap(permK);
	_RHSPermutation[domain].swap(permRHS);

	ROW.push_back(COL.size() + 1);

	_data->K[domain].resize(_domainDOFsSize[domain], _domainDOFsSize[domain], COL.size());
	_data->K[domain].fillPattern(_domainDOFsSize[domain], ROW.data(), COL.data());
	_data->K[domain].structureUpdated();
}
void NodesUniformFETIComposer::_buildKBEMPattern(esint domain)
{
//	SurfaceStore* surface = info::mesh->domainsSurface;
//	_data->K[domain].rows = surface->cdistribution[domain + 1] - surface->cdistribution[domain];
//	_data->K[domain].cols = surface->cdistribution[domain + 1] - surface->cdistribution[domain];
//	_data->K[domain].mtype = MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
//	_data->K[domain].type = 'S';
//	_data->f[domain].resize(data->K[domain].rows);
//
//	_data->K[domain].CSR_I_row_indices.push_back(1);
//	for (esint r = 0; r < data->K[domain].rows; r++) {
//		for (esint c = r; c < data->K[domain].cols; c++) {
//			data->K[domain].CSR_J_col_indices.push_back(c + 1);
//		}
//		data->K[domain].CSR_I_row_indices.push_back(data->K[domain].CSR_J_col_indices.size() + 1);
//	}
//
//	_data->K[domain].CSR_V_values.resize(data->K[domain].CSR_J_col_indices.size());
//	_data->K[domain].nnz = data->K[domain].CSR_J_col_indices.size();
//
//	esint RHSsize = 0;
//	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
//		if (
//				info::mesh->boundaryRegions[r]->dimension &&
//				info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
//
//			esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
//			esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
//			auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;
//
//			for (esint i = begin; i < end; ++i, ++enodes) {
//				RHSsize += enodes->size() * _DOFs;
//			}
//		}
//	}
//
//	std::vector<esint> permRHS(RHSsize);
//	std::vector<esint> RHSPattern(RHSsize), ROW, COL;
//
//	esint *RHSoffset = RHSPattern.data();
//
//	auto insert = [&] (serializededata<esint, esint>::const_iterator &enodes) {
//		esint *_RHS = RHSoffset;
//		for (auto n = enodes->begin(); n != enodes->end(); ++n, ++RHSoffset) {
//			auto DOFs = (_DOFMap->begin() + *n)->begin();
//			while (*DOFs != domain + info::mesh->elements->firstDomain) {
//				DOFs += 1 + _DOFs;
//			}
//			*RHSoffset = *(DOFs + 1);
//		}
//		for (int dof = 1; dof < _DOFs; ++dof) {
//			for (size_t n = 0; n < enodes->size(); ++n, ++RHSoffset) {
//				*RHSoffset = *(_RHS + n) + dof;
//			}
//		}
//	};
//
//	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
//		if (
//				info::mesh->boundaryRegions[r]->dimension &&
//				info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
//
//			esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
//			esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
//			auto enodes = info::mesh->boundaryRegions[r]->procNodes->cbegin() + begin;
//
//			for (esint i = begin; i < end; ++i, ++enodes) {
//				insert(enodes);
//			}
//		}
//	}
//
//	std::vector<esint> pRHS(RHSPattern.size());
//	std::iota(pRHS.begin(), pRHS.end(), 0);
//	std::sort(pRHS.begin(), pRHS.end(), [&] (esint i, esint j) {
//		return RHSPattern[i] < RHSPattern[j];
//	});
//
//	if (pRHS.size()) {
//		permRHS[pRHS.front()] = 0;
//		for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
//			if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
//				++nonzeros;
//			}
//			permRHS[pRHS[i]] = nonzeros;
//		}
//	}
//
//	_RHSPermutation[domain].swap(permRHS);
}

void NodesUniformFETIComposer::_buildDirichlet()
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

void NodesUniformFETIComposer::_buildMortars()
{
	struct __lambda__ { esint lambda, pair, id; };

	esint lambdaCount = 0;
	std::vector<__lambda__> lmap;
	auto &myids = info::mesh->nodes->IDs->datatarray();

	// collect my nodes in B in order to compute mortar lambdas
	if (info::mesh->contacts->B.size()) {
		for (auto it = info::mesh->contacts->B.begin(), prev = it; it != info::mesh->contacts->B.end(); ++it) {
			if (it == info::mesh->contacts->B.begin() || prev->pair != it->pair || prev->i != it->i) {
				auto nit = std::lower_bound(myids.begin() + info::mesh->nodes->uniqInfo.nhalo, myids.end(), it->i);
				if (nit != myids.end() && *nit == it->i) {
					lmap.push_back({ lambdaCount++, it->pair, it->i });
				}
				prev = it;
			}
		}
	}

	esint size = Communication::exscan(lambdaCount);
	for (auto lit = lmap.begin(); lit != lmap.end(); ++lit) {
		lit->lambda += lambdaCount;
	}

	std::vector<std::vector<__lambda__> > sLambdas(info::mesh->neighborsWithMe.size(), lmap), rLambdas(info::mesh->neighborsWithMe.size());
	if (!Communication::exchangeUnknownSize(sLambdas, rLambdas, info::mesh->neighborsWithMe)) {
		eslog::internalFailure("cannot exchange mortar lambdas.\n");
	}

	// lambdas <node, lambda> are sorted according to lambda in order to get correct result from FETI solver
	std::vector<__lambda__> lambdas;
	for (size_t n = 0; n < rLambdas.size(); ++n) {
		lambdas.insert(lambdas.end(), rLambdas[n].begin(), rLambdas[n].end());
	}
	std::sort(lambdas.begin(), lambdas.end(), [] (const __lambda__ &l1, const __lambda__ &l2) { return l1.lambda < l2.lambda; });

//	Communication::serialize([&] () {
//		std::cout << info::mpi::rank << ": \n";
//		std::cout << lambdas;
////		auto dof = _DOFMap->begin();
////		for (esint n = 0; n < info::mesh->nodes->size; ++n, ++dof) {
////			std::cout << info::mesh->nodes->IDs->datatarray()[n] << "[" << info::mesh->nodes->coordinates->datatarray()[n] << "]: " << *dof << "\n";
////		}
//	});

	// found nodes that are on mortar interface -> we have exchange domain map in order to get correct decomposition
	std::vector<esint> mids, mylambdas;
	for (auto it = info::mesh->contacts->B.begin(); it != info::mesh->contacts->B.end(); ++it) {
		auto nit = std::find(myids.begin(), myids.begin() + info::mesh->nodes->uniqInfo.nhalo, it->j);
		if (nit == myids.begin() + info::mesh->nodes->uniqInfo.nhalo || *nit != it->j) {
			nit = std::lower_bound(myids.begin() + info::mesh->nodes->uniqInfo.nhalo, myids.end(), it->j);
		}
		if (nit != myids.end() && *nit == it->j) {
			mids.push_back(it->j);
			if (mylambdas.size() == 0 || mylambdas.back() != it->i) {
				mylambdas.push_back(it->i);
			}
		}
	}
	utils::sortAndRemoveDuplicates(mylambdas);
	utils::sortAndRemoveDuplicates(mids);
//	Communication::serialize([&] () {
//		std::cout << info::mpi::rank << ": \n";
//		std::cout << mylambdas;
//	});

	std::vector<esint> sBuffer;
	std::vector<std::vector<esint> > rBuffer(info::mesh->neighborsWithMe.size());
	sBuffer.push_back(mids.size());
	for (auto it = mids.begin(); it != mids.end(); ++it) {
		auto nit = std::find(myids.begin(), myids.begin() + info::mesh->nodes->uniqInfo.nhalo, *it);
		if (nit == myids.begin() + info::mesh->nodes->uniqInfo.nhalo || *nit != *it) {
			nit = std::lower_bound(myids.begin() + info::mesh->nodes->uniqInfo.nhalo, myids.end(), *it);
		}
		if (nit != myids.end() && *nit == *it) {
			sBuffer.push_back(*it);
			size_t size = sBuffer.size();
			sBuffer.push_back(0);
			auto dmap = _DOFMap->begin() + (nit - myids.begin());
			for (auto di = dmap->begin(); di != dmap->end(); ++di) {
				if (_data->K.ismy(di->domain)) {
					++sBuffer[size];
					sBuffer.push_back(di->domain);
					sBuffer.push_back(di->index);
				}
			}
		}
	}

	if (!Communication::exchangeUnknownSize(sBuffer, rBuffer, info::mesh->neighborsWithMe)) {
		eslog::internalFailure("cannot exchange mortar d-map.\n");
	}

	struct npair { esint id, n, offset; };
	std::vector<npair> moffset;
	for (size_t n = 0; n < rBuffer.size(); ++n) {
		for (size_t i = 1; i < rBuffer[n].size(); i += 2 * rBuffer[n][i + 1] + 2) {
			moffset.push_back(npair{ rBuffer[n][i], (esint)n, (esint)i });
		}
	}
	std::sort(moffset.begin(), moffset.end(), [&] (const npair &i, const npair &j) {
		if (i.id == j.id) {
			return i.n < j.n;
		}
		return i.id < j.id;
	});

	std::vector<esint> didist = { 0 };
	std::vector<DI> didata;

	std::vector<std::vector<esint> > rows(info::mesh->elements->ndomains), cols(info::mesh->elements->ndomains);
	std::vector<std::vector<double> > vals(info::mesh->elements->ndomains);

	// build B1 in correct order
	for (auto lambda = lambdas.begin(); lambda != lambdas.end(); ++lambda) {
		if (!std::binary_search(mylambdas.begin(), mylambdas.end(), lambda->id)) {
			continue;
		}
		auto begin = std::lower_bound(info::mesh->contacts->B.begin(), info::mesh->contacts->B.end(), *lambda, [] (const ijv &b, const __lambda__ &l) { return b.pair == l.pair ? b.i < l.id : b.pair < l.pair; });
		auto end = begin;
		while (end != info::mesh->contacts->B.end() && lambda->pair == end->pair && lambda->id == end->i) {
			++end;
		}

		for (auto it = begin; it != end; ++it) {
			auto nit = std::lower_bound(moffset.begin(), moffset.end(), it->j, [&] (const npair &info, const esint &lambda) { return info.id < lambda; });
			esint domains = 0;
			for (auto nnit = nit; nnit != moffset.end() && nnit->id == it->j; ++nnit) {
				domains += rBuffer[nnit->n][nnit->offset + 1];
			}
			while (nit != moffset.end() && nit->id == it->j) {
				if (info::mesh->neighborsWithMe[nit->n] == info::mpi::rank) {
					for (esint i = 0; i < rBuffer[nit->n][nit->offset + 1]; ++i) {
						rows[rBuffer[nit->n][nit->offset + 2 + 2 * i] - _data->K.doffset].push_back(lambda->lambda + 1);
						cols[rBuffer[nit->n][nit->offset + 2 + 2 * i] - _data->K.doffset].push_back(rBuffer[nit->n][nit->offset + 2 + 2 * i + 1] + 1);
						vals[rBuffer[nit->n][nit->offset + 2 + 2 * i] - _data->K.doffset].push_back(it->v / domains);
					}
				}
				for (esint i = 0; i < rBuffer[nit->n][nit->offset + 1]; ++i) {
					didata.push_back(DI{ rBuffer[nit->n][nit->offset + 2 + 2 * i], rBuffer[nit->n][nit->offset + 2 + 2 * i + 1] });
				}
				++nit;
			}
		}
		if (begin != end) {
			std::sort(didata.begin() + didist.back(), didata.end());
			didist.push_back(didata.size());
		}
	}

	serializededata<esint, DI> *dmap = new serializededata<esint, DI>(tarray<esint>(info::env::OMP_NUM_THREADS, didist), tarray<DI>(info::env::OMP_NUM_THREADS, didata));

//	Communication::serialize([&] () {
//		std::cout << info::mpi::rank << "<" << info::mesh->elements->firstDomain << "> " << ": \n";
//		for (size_t n = 0; n < rows.size(); ++n) {
//			for (size_t i = 0; i < rows[n].size(); i++) {
//				if (i && rows[n][i] != rows[n][i - 1]) {
//					printf("\n");
//				}
//				printf("%lu:%2d: <%lu:%d> = %8.6f\n", n, rows[n][i], info::mesh->elements->firstDomain + n, cols[n][i] - 1, vals[n][i]);
//			}
//			printf("\n");
//		}
////		for (esint n = 0; n < info::mesh->nodes->size; ++n) {
////			std::cout << "n[" << info::mesh->nodes->IDs->datatarray()[n] << "]: " << *(_DOFMap->begin() + n) << "\n";
////		}
////		std::cout << *_DOFMap << "\n\n";
//		std::cout << lambdas.size() << " vs. " << dmap->structures() << "\n";
////		auto ddit = dmap->begin();
////		for (auto ll = lambdas.begin(); ll != lambdas.end(); ++ll, ++ddit) {
////			std::cout << *ll << ": " << *ddit << "\n";
////		}
////		std::cout << *dmap << "\n\n";
////		std::cout << *info::mesh->nodes->coordinates << "\n";
////		std::cout << *_DOFMap << "\n";
////		for (esint n = 0; n < info::mesh->nodes->size; ++n) {
////			std::cout << info::mesh->nodes->IDs->datatarray()[n] << "::" << info::mesh->nodes->coordinates->datatarray()[n] << " -> " << *(_DOFMap->begin() + n) << "\n";
////		}
//	});



	_data->mortars.initDomains(info::mesh->elements->ndomains);
	_data->mortars.fillDecomposition(
				info::mpi::rank, info::mpi::size, info::mesh->neighbors.size(),
				_data->K.distribution, info::mesh->neighbors.data(), dmap);

	for (esint d = 0; d < info::mesh->elements->ndomains; ++d) {
		_data->mortars[d].resize(size, _data->K[d].ncols, vals[d].size());
		_data->mortars[d].fillPattern(rows[d].size(), rows[d].data(), cols[d].data());
		_data->mortars[d].fillValues(vals[d].size(), vals[d].data());
	}
	_data->mortars.structureUpdated();

	delete dmap;
}

void NodesUniformFETIComposer::_buildInequality()
{
	std::vector<std::pair<esint, esint> > indices;
	kernel->solverDataProvider->general->inequalityIndices(indices);

	std::vector<esint> dofs, permutation, pattern;
	for (auto i = indices.begin(); i != indices.end(); ++i) {
		dofs.push_back(i->first * _DOFs + i->second);
	}
	permutation.resize(dofs.size());
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&] (esint i, esint j) {
		return dofs[i] < dofs[j];
	});
	std::sort(dofs.begin(), dofs.end());
	pattern = dofs;
	utils::removeDuplicates(pattern);

	_data->gapDirection.initVectors(1);
	_data->gapDirection.resize(info::mesh->nodes->size * _DOFs, indices.size());
	_data->gapDirection.holder()->fillPattern(pattern.data());
	_data->gap.shallowCopyStructure(&_data->gapDirection);

	std::vector<double> values(indices.size());
	kernel->solverDataProvider->general->inequalityNormals(values);
	fillPermutedSparseData(_data->gapDirection[0].vals, dofs, permutation, values);

	kernel->solverDataProvider->general->inequalityGaps(values);
	fillPermutedSparseData(_data->gap[0].vals, dofs, permutation, values);
}

void NodesUniformFETIComposer::synchronize(const Builder &builder)
{

}

void NodesUniformFETIComposer::computeCornerNodes()
{
	if (info::mesh->nodes->domains == NULL) {
		mesh::computeNodeDomainDistribution();
	}

	std::vector<esint> uniq;
	std::vector<std::vector<esint> > nodes;

	esint i = 0;
	for (auto dmap = info::mesh->nodes->domains->cbegin(); dmap != info::mesh->nodes->domains->cend(); ++dmap, ++i) {
		if (dmap->size() > 1) {
			size_t index = 0;
			for (auto u = uniq.begin(); u != uniq.end(); ++u, ++index) {
				if (*u == (esint)dmap->size() && memcmp(&*(u + 1), dmap->data(), dmap->size() * sizeof(esint)) == 0) {
					break;
				}
				u += *u;
			}
			if (index == nodes.size()) {
				uniq.push_back(dmap->size());
				uniq.insert(uniq.end(), dmap->begin(), dmap->end());
				nodes.push_back({ i });
			} else {
				nodes[index].push_back(i);
			}
		}
	}

	for (size_t n = 0; n < nodes.size(); ++n) {
		esint inc = nodes[n].size() / 3;
		inc = inc ? inc : 1;
		for (size_t nn = 0; nn < nodes[n].size(); nn += inc) {
			info::mesh->FETIData->corners.push_back(nodes[n][nn]);
		}
	}

	std::sort(info::mesh->FETIData->corners.begin(), info::mesh->FETIData->corners.end());

	eslog::checkpointln("MESH: CORNER NODES COMPUTED");
}

static void addFixPoints(const serializededata<esint, esint>* elements, esint begin, esint end, const serializededata<esint, Element*>* epointers, std::vector<esint> &fixPoints)
{
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
	auto element = elements->begin() + begin;
	const auto &epointer = epointers->datatarray();
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
	for (size_t i = 0; i < permutation.size(); i++) {
		if (i && originnodes[permutation[i]] != originnodes[permutation[i - 1]]) {
			utils::sortAndRemoveDuplicates(data, dist.back());
			dist.push_back(data.size());
			ids.push_back(originnodes[permutation[i]]);
		}
		data.push_back(neighsnodes[permutation[i]]);
	}
	utils::sortAndRemoveDuplicates(data, dist.back());
	dist.push_back(data.size());

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
			float last_l = pids[p].size(), l = 1;

			while (fabs((l - last_l) / l) > 1e-6) {
				MATH::upCSRMatVecProduct(pids[p].size(), pids[p].size(), pdist[p].data(), pdata[p].data(), vals.data(), x.data(), y.data());
				last_l = l;
				l = MATH::vecNorm(pids[p].size(), y.data());
				MATH::vecScale(pids[p].size(), 1 / l, y.data());
				x.swap(y);
			}

			fixPoints.push_back(ids[pids[p][MATH::vecNormMaxIndex(pids[p].size(), x.data())]]);
		}
	}
}

void NodesUniformFETIComposer::computeFixPoints()
{
	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > fixPoints(threads), fixPointsDist(threads);

	fixPointsDist.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> dist, data, partition;
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			size_t size = fixPoints[t].size();
			addFixPoints(info::mesh->elements->procNodes, info::mesh->elements->elementsDistribution[d], info::mesh->elements->elementsDistribution[d + 1], info::mesh->elements->epointers, fixPoints[t]);
			utils::sortAndRemoveDuplicates(fixPoints[t], size);
			fixPointsDist[t].push_back(fixPoints[t].size());
		}
	}

	utils::threadDistributionToFullDistribution(fixPointsDist);
	info::mesh->FETIData->iFixPointsDistribution.clear();
	for (size_t t = 0; t < threads; t++) {
		info::mesh->FETIData->iFixPointsDistribution.insert(info::mesh->FETIData->iFixPointsDistribution.end(), fixPointsDist[t].begin(), fixPointsDist[t].end());
		info::mesh->FETIData->innerFixPoints.insert(info::mesh->FETIData->innerFixPoints.end(), fixPoints[t].begin(), fixPoints[t].end());
	}

	DebugOutput::surfaceFixPoints();
}

void computeFixPointsOnSurface()
{
	if (info::mesh->domainsSurface == NULL || info::mesh->domainsSurface->enodes == NULL) {
		mesh::computeDomainsSurface();
	}

	size_t threads = info::env::OMP_NUM_THREADS;

	std::vector<std::vector<esint> > fixPoints(threads), fixPointsDist(threads);

	fixPointsDist.front().push_back(0);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<esint> dist, data, partition;
		for (esint d = info::mesh->elements->domainDistribution[t]; d < info::mesh->elements->domainDistribution[t + 1]; d++) {
			size_t size = fixPoints[t].size();
			addFixPoints(info::mesh->domainsSurface->enodes, info::mesh->domainsSurface->edistribution[d], info::mesh->domainsSurface->edistribution[d + 1], info::mesh->domainsSurface->epointers, fixPoints[t]);
			utils::sortAndRemoveDuplicates(fixPoints[t], size);
			fixPointsDist[t].push_back(fixPoints[t].size());
		}
	}

	utils::threadDistributionToFullDistribution(fixPointsDist);
	info::mesh->FETIData->sFixPointsDistribution.clear();
	for (size_t t = 0; t < threads; t++) {
		info::mesh->FETIData->sFixPointsDistribution.insert(info::mesh->FETIData->sFixPointsDistribution.end(), fixPointsDist[t].begin(), fixPointsDist[t].end());
		info::mesh->FETIData->surfaceFixPoints.insert(info::mesh->FETIData->surfaceFixPoints.end(), fixPoints[t].begin(), fixPoints[t].end());
	}

	DebugOutput::surfaceFixPoints();
}

void NodesUniformFETIComposer::buildB0FromCorners(MatrixIJVFETI &B0)
{
	if (info::mesh->FETIData->corners.size() == 0) {

	}

	esint lambda = 0;
	std::vector<std::vector<esint> > ROWS(info::mesh->elements->ndomains), COLS(info::mesh->elements->ndomains);
	std::vector<std::vector<double> > VALS(info::mesh->elements->ndomains);

	auto dmap = _DOFMap->begin();
	for (size_t n = 0, prev = 0; n < info::mesh->FETIData->corners.size(); prev = n++) {
		dmap += _DOFs * (info::mesh->FETIData->corners[n] - info::mesh->FETIData->corners[prev]);
		for (int dof = 0; dof < _DOFs; ++dof) {
			for (auto di1 = (dmap + dof)->begin(), di2 = di1 + 1; di2 != (dmap + dof)->end(); ++di1, ++di2) {
				if (_data->K.ismy(di1->domain) && _data->K.ismy(di2->domain)) {
					ROWS[di1->domain - _data->K.doffset].push_back(lambda + 1);
					COLS[di1->domain - _data->K.doffset].push_back(di1->index + 1);
					VALS[di1->domain - _data->K.doffset].push_back(1);

					ROWS[di2->domain - _data->K.doffset].push_back(lambda + 1);
					COLS[di2->domain - _data->K.doffset].push_back(di2->index + 1);
					VALS[di2->domain - _data->K.doffset].push_back(-1);
					++lambda;
				}
			}
		}
	}

	_data->B0.initDomains(_data->K.domains);
	#pragma omp parallel for
	for (esint d = 0; d < info::mesh->elements->ndomains; d++) {
		_data->B0[d].type = MatrixType::REAL_UNSYMMETRIC;
		_data->B0[d].resize(lambda, _data->K[d].ncols, VALS[d].size());
		_data->B0[d].fillPattern(ROWS[d].size(), ROWS[d].data(), COLS[d].data());
		_data->B0[d].fillValues(VALS[d].size(), VALS[d].data());
	}
}




