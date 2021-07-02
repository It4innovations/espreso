
#include "nodes.uniform.feti.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/nodestore.h"

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

UniformNodesFETIPattern::UniformNodesFETIPattern(): dmap(nullptr)
{

}

UniformNodesFETIPattern::~UniformNodesFETIPattern()
{
	if (dmap) { delete dmap; }
}

void fillDomainMap(UniformNodesFETIPattern *pattern, int dofs)
{
	pattern->elements.resize(info::mesh->domains->size);

	// nID, domain
	std::vector<std::vector<std::pair<esint, esint> > > ntodomains(info::env::threads);

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; t++) {
		std::vector<std::pair<esint, esint> > tdata;

		for (size_t d = info::mesh->domains->distribution[t]; d != info::mesh->domains->distribution[t + 1]; ++d) {
			std::vector<esint> dnodes(
					(info::mesh->elements->nodes->begin() + info::mesh->domains->elements[d])->begin(),
					(info::mesh->elements->nodes->begin() + info::mesh->domains->elements[d + 1])->begin());

			utils::sortAndRemoveDuplicates(dnodes);
			for (size_t i = 0; i < dnodes.size(); i++) {
				tdata.push_back(std::pair<esint, esint>(dnodes[i], d));
			}
			pattern->elements[d].size = dnodes.size() * dofs;
		}

		ntodomains[t].swap(tdata);
	}

	for (int t = 1; t < info::env::threads; t++) {
		ntodomains[0].insert(ntodomains[0].end(), ntodomains[t].begin(), ntodomains[t].end());
	}

	std::sort(ntodomains[0].begin(), ntodomains[0].end());

	std::vector<std::vector<esint> > DOFs(info::env::threads);
	std::vector<size_t> ndistribution = tarray<size_t>::distribute(info::env::threads, ntodomains[0].size());

	for (int t = 1; t < info::env::threads; t++) {
		while (
				ndistribution[t] < ntodomains[0].size() &&
				ntodomains[0][ndistribution[t]].first == ntodomains[0][ndistribution[t] - 1].first) {

			++ndistribution[t];
		}
	}

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; t++) {
		std::vector<esint> tDOFs(info::mesh->domains->size);

		for (size_t n = ndistribution[t]; n < ndistribution[t + 1]; ++n) {
			tDOFs[ntodomains[0][n].second] += dofs;
		}

		DOFs[t].swap(tDOFs);
	}

	utils::sizesToOffsets(DOFs);

	std::vector<std::vector<std::vector<esint> > > sBuffer(info::env::threads);
	std::vector<std::vector<esint> > rBuffer(info::mesh->neighborsWithMe.size());

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; t++) {
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
					DOFs[t][ntodomains[0][i].second] += dofs;
				}
			}
		}
		sBuffer[t].swap(tBuffer);
	}

	for (int t = 1; t < info::env::threads; t++) {
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

		for (int dof = 0; dof < dofs; ++dof) {
			for (size_t i = 0; i < current.size(); ++i) {
				DOFData.push_back({ current[i].domain, current[i].index + dof });
			}
			DOFDistribution.push_back(DOFData.size());
		}
		current.clear();
	}

	std::vector<size_t> distribution = info::mesh->nodes->distribution, datadistribution(info::env::threads + 1);
	for (int t = 1; t < info::env::threads; t++) {
		distribution[t] = dofs * distribution[t] + 1;
		datadistribution[t] = distribution[t] < DOFDistribution.size() ? DOFDistribution[distribution[t]] : DOFDistribution.back();
	}
	datadistribution[info::env::threads] = dofs * distribution[info::env::threads] < DOFDistribution.size() ? DOFDistribution[dofs * distribution[info::env::threads]] : DOFDistribution.back();
	distribution[info::env::threads] = dofs * distribution[info::env::threads] + 1;

	pattern->dmap = new serializededata<esint, DI>(
			tarray<esint>(distribution, DOFDistribution),
			tarray<DI>(datadistribution, DOFData));
}

void fillPermutation(UniformNodesFETIPattern *pattern, int dofs, Matrix_Type type, int domain)
{
	auto sizeUpper = [] (int size) {
		return (size * size - size) / 2 + size;
	};

	auto sizeFull = [] (int size) {
		return size * size;
	};

	auto fillUpper = [] (IJ* &target, const esint *begin, const esint *end) {
		for (auto row = begin, colbegin = begin; row != end; ++row, ++colbegin) {
			for (auto col = colbegin; col != end; ++col, ++target) {
				if (*row <= *col) {
					target->row = *row;
					target->column = *col;
				} else {
					target->row = *col;
					target->column = *row;
				}
			}
		}
	};

	auto fillFull = [] (IJ* &target, const esint *begin, const esint *end) {
		for (auto row = begin; row != end; ++row) {
			for (auto col = begin; col != end; ++col, ++target) {
				target->row = *row;
				target->column = *col;
			}
		}
	};

	std::function<int(int)> size;
	std::function<void(IJ*&, const esint*, const esint*)> fill;

	switch (type) {
	case Matrix_Type::REAL_SYMMETRIC_INDEFINITE: size = sizeUpper; fill = fillUpper; break;
	case Matrix_Type::REAL_SYMMETRIC_POSITIVE_DEFINITE: size = sizeUpper; fill = fillUpper; break;
	case Matrix_Type::REAL_UNSYMMETRIC: size = sizeFull; fill = fillFull; break;
	}

	auto ebegin = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain];
	auto eend = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain + 1];

	esint Ksize = 0, RHSsize = 0;
	for (auto e = ebegin; e != eend; ++e) {
		RHSsize += e->size() * dofs;
		Ksize += size(e->size() * dofs);
	}

	std::vector<IJ> KPattern(Ksize);
	std::vector<esint> RHSPattern(RHSsize);

	IJ *Koffset = KPattern.data();
	esint *RHSoffset = RHSPattern.data();

	for (auto e = ebegin; e != eend; ++e) {
		esint *_RHS = RHSoffset;
		for (int dof = 0; dof < dofs; ++dof) {
			for (auto n = e->begin(); n != e->end(); ++n, ++RHSoffset) {
				auto DOFs = (pattern->dmap->begin() + (*n * dofs + dof))->begin();
				while (DOFs->domain != domain + info::mesh->domains->offset) {
					++DOFs;
				}
				*RHSoffset = DOFs->index;
			}
		}
		fill(Koffset, _RHS, RHSoffset);
	}

	std::vector<IJ> KData = KPattern;
	std::vector<esint> RHSData = RHSPattern;

	utils::sortAndRemoveDuplicates(KPattern);
	utils::sortAndRemoveDuplicates(RHSPattern);

	pattern->elements[domain].row.reserve(KPattern.size());
	pattern->elements[domain].column.reserve(KPattern.size());
	pattern->elements[domain].K.reserve(KData.size());
	pattern->elements[domain].f.reserve(RHSData.size());

	for (size_t i = 0; i < KPattern.size(); ++i) {
		pattern->elements[domain].row.push_back(KPattern[i].row);
		pattern->elements[domain].column.push_back(KPattern[i].column);
	}

	for (size_t i = 0; i < KData.size(); ++i) {
		pattern->elements[domain].K.push_back(std::lower_bound(KPattern.begin(), KPattern.end(), KData[i]) - KPattern.begin());
	}
	for (size_t i = 0; i < RHSData.size(); ++i) {
		pattern->elements[domain].f.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), RHSData[i]) - RHSPattern.begin());
	}
}

void UniformNodesFETIPattern::set(int dofs, Matrix_Type type)
{
	if (dmap) { delete dmap; }

	fillDomainMap(this, dofs);
	for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
		fillPermutation(this, dofs, type, domain);
	}
}

void UniformNodesFETIPattern::fillCSR(int domain, esint *rows, esint *cols)
{
	rows[0] = cols[0] = _Matrix_CSR_Pattern::Indexing;
	for (size_t r = 1, c = 1; c < elements[domain].column.size(); ++c) {
		cols[c] = elements[domain].column[c] + _Matrix_CSR_Pattern::Indexing;
		if (elements[domain].row[c] != elements[domain].row[c - 1]) {
			rows[r++] = c + _Matrix_CSR_Pattern::Indexing;
		}
	}
}
