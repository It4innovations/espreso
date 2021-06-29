
#include "nodes.uniform.feti.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "wrappers/mpi/communication.h"

#include <numeric>

using namespace espreso;

struct IJ {
	esint row, column;
};

inline bool operator==(const IJ &left, const IJ &right)
{
	return left.row == right.row && left.column == right.column;
}

inline bool operator!=(const IJ &left, const IJ &right)
{
	return !(left == right);
}

inline bool operator<(const IJ &left, const IJ &right)
{
	return left.row == right.row ? left.column < right.column : left.row < right.row;
}


void UniformNodesFETIPattern::init()
{
	// TODO: improve performance using info::mesh->nodes->domains
	size_t threads = info::env::OMP_NUM_THREADS;

	elements.resize(info::mesh->domains->size);

	// nID, domain
	std::vector<std::vector<std::pair<esint, esint> > > ntodomains(threads);

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<std::pair<esint, esint> > tdata;

		for (size_t d = info::mesh->domains->distribution[t]; d != info::mesh->domains->distribution[t + 1]; ++d) {
			std::vector<esint> dnodes(
					(info::mesh->elements->nodes->begin() + info::mesh->domains->elements[d])->begin(),
					(info::mesh->elements->nodes->begin() + info::mesh->domains->elements[d + 1])->begin());

			utils::sortAndRemoveDuplicates(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				tdata.push_back(std::pair<esint, esint>(dnodes[i], d));
			}
			elements[d].size = dnodes.size() * _DOFs;
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
		std::vector<esint> tDOFs(info::mesh->domains->size);

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
		std::vector<std::vector<esint> > tBuffer(info::mesh->neighborsWithMe.size());

		if (ndistribution[t] < ndistribution[t  +1]) {
			auto nranks = info::mesh->nodes->ranks->begin() + ntodomains[0][ndistribution[t]].first;
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
						tBuffer[noffset].push_back(info::mesh->domains->offset + ntodomains[0][i].second);
						tBuffer[noffset].push_back(DOFs[t][ntodomains[0][i].second]);
					}
				}
				++nranks;

				for (size_t i = begin; i < n; i++) {
					DOFs[t][ntodomains[0][i].second] += _DOFs;
				}
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
		datadistribution[t] = distribution[t] < DOFDistribution.size() ? DOFDistribution[distribution[t]] : DOFDistribution.back();
	}
	datadistribution[threads] = _DOFs * distribution[threads] < DOFDistribution.size() ? DOFDistribution[_DOFs * distribution[threads]] : DOFDistribution.back();
	distribution[threads] = _DOFs * distribution[threads] + 1;

	_DOFMap = new serializededata<esint, DI>(
			tarray<esint>(distribution, DOFDistribution),
			tarray<DI>(datadistribution, DOFData));
}

void init(esint domain, int dofs, serializededata<esint, DI> *map, std::function<int(int)> getMatrixSize, std::function<void(IJ *target, const esint *begin, const esint *end)> insertPatter)
{
	auto ebegin = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain];
	auto eend = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain + 1];

	esint Ksize = 0, RHSsize = 0;
	for (auto e = ebegin; e != eend; ++e) {
		RHSsize += e->size() * dofs;
		Ksize += getMatrixSize(e->size() * dofs);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			if (info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
				esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
				esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
				auto enodes = info::mesh->boundaryRegions[r]->elements->cbegin() + begin;

				for (esint i = begin; i < end; ++i, ++enodes) {
					RHSsize += enodes->size() * dofs;
					Ksize += getMatrixSize(enodes->size() * dofs);
				}
			}
		}
	}
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (!info::mesh->boundaryRegions[r]->dimension) {
			esint prev = 0;
			auto dmap = map->begin();
			for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); prev = *n++) {
				dmap += (*n - prev) * dofs;
				if (dmap->begin()->domain == domain + info::mesh->domains->offset) {
					RHSsize += dofs;
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
		for (int dof = 0; dof < dofs; ++dof) {
			for (auto n = enodes->begin(); n != enodes->end(); ++n, ++RHSoffset) {
				auto DOFs = (map->begin() + (*n * dofs + dof))->begin();
				while (DOFs->domain != domain + info::mesh->domains->offset) {
					++DOFs;
				}
				*RHSoffset = DOFs->index;
			}
		}
		insertPatter(Koffset, _RHS, RHSoffset);
	};

	for (auto e = ebegin; e != eend; ++e) {
		insert(e);
		Koffset += getMatrixSize(e->size() * dofs);
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			if (info::mesh->boundaryRegions[r]->eintervalsDistribution[domain] < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]) {
				esint begin = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]].begin;
				esint end = info::mesh->boundaryRegions[r]->eintervals[info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1] - 1].end;
				auto enodes = info::mesh->boundaryRegions[r]->elements->cbegin() + begin;
				for (esint i = begin; i < end; ++i, ++enodes) {
					insert(enodes);
					Koffset += getMatrixSize(enodes->size() * dofs);
				}
			}
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); r++) {
		if (!info::mesh->boundaryRegions[r]->dimension) {
			esint prev = 0;
			auto dmap = map->begin();
			for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().begin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); prev = *n++) {
				dmap += (*n - prev) * dofs;
				for (int dof = 0; dof < dofs; ++dof) {
					if ((dmap + dof)->begin()->domain == domain + info::mesh->domains->offset) {
						*RHSoffset++ = (dmap + dof)->begin()->index;
					}
				}
			}
		}
	}

//	std::vector<esint> pK(KPattern.size());
//	std::iota(pK.begin(), pK.end(), 0);
//	std::sort(pK.begin(), pK.end(), [&] (esint i, esint j) {
//		return KPattern[i] < KPattern[j];
//	});
//
//	std::vector<esint> pRHS(RHSPattern.size());
//	std::iota(pRHS.begin(), pRHS.end(), 0);
//	std::sort(pRHS.begin(), pRHS.end(), [&] (esint i, esint j) {
//		return RHSPattern[i] < RHSPattern[j];
//	});
//
//	ROW.reserve(_domainDOFsSize[domain] + 1);
//	COL.reserve(KPattern.size());
//	ROW.push_back(1);
//	COL.push_back(KPattern[pK.front()].column + 1);
//	permK[pK.front()] = 0;
//	for (size_t i = 1, nonzeros = 0; i < pK.size(); i++) {
//		if (KPattern[pK[i]] != KPattern[pK[i - 1]]) {
//			++nonzeros;
//			COL.push_back(KPattern[pK[i]].column + 1);
//			if (KPattern[pK[i - 1]].row != KPattern[pK[i]].row) {
//				ROW.push_back(nonzeros + 1);
//			}
//		}
//		permK[pK[i]] = nonzeros;
//	}
//	permRHS[pRHS.front()] = 0;
//	for (size_t i = 1, nonzeros = 0; i < pRHS.size(); i++) {
//		if (RHSPattern[pRHS[i]] != RHSPattern[pRHS[i - 1]]) {
//			++nonzeros;
//		}
//		permRHS[pRHS[i]] = nonzeros;
//	}
//
//	_KPermutation[domain].swap(permK);
//	_RHSPermutation[domain].swap(permRHS);
//
//	ROW.push_back(COL.size() + 1);
//
//	_data->K[domain].resize(_domainDOFsSize[domain], _domainDOFsSize[domain], COL.size());
//	_data->K[domain].fillPattern(_domainDOFsSize[domain], ROW.data(), COL.data());
//	_data->K[domain].structureUpdated();
}

void UniformNodesFETIPattern::initUpper()
{
//	for (esint d = 0; d < info::mesh->domains->size; ++d) {
//		init(d, [] (int dofs) { return (dofs * dofs - dofs) / 2 + dofs; });
//	}
}

void UniformNodesFETIPattern::initFull()
{
//	for (esint d = 0; d < info::mesh->domains->size; ++d) {
//		init(d, [] (int dofs) { return dofs * dofs; });
//	}
}
