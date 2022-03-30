
#include "nodes.uniform.feti.h"
#include "ij.h"

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

#include "basis/utilities/print.h"

using namespace espreso;

UniformNodesFETIPattern::UniformNodesFETIPattern()
: dofs(0)
{

}

UniformNodesFETIPattern::~UniformNodesFETIPattern()
{

}

void fillDecomposition(UniformNodesFETIPattern *pattern, int dofs, DOFsDecomposition &decomposition)
{
	pattern->elements.resize(info::mesh->domains->size);
	pattern->bregion.resize(info::mesh->domains->size);

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
	std::vector<DIndex> DOFData;

	// TODO: make it parallel
	// parallelization is possible if node order will be kept as: boundary first!
	// now we prefer generality
	auto nranks = info::mesh->nodes->ranks->begin();
	std::vector<esint> roffset(rBuffer.size());
	std::vector<DIndex> current; current.reserve(50);
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

	decomposition.dbegin = info::mesh->domains->offset;
	decomposition.dend = info::mesh->domains->offset + info::mesh->domains->size;
	decomposition.dtotal = info::mesh->domains->totalSize;
	decomposition.neighDomain.resize(info::mesh->neighbors.size() + 1);
	auto ddist = info::mesh->domains->gatherProcDistribution(); // remove this
	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		decomposition.neighDomain[n] = ddist[info::mesh->neighbors[n]];
	}
	decomposition.neighDomain.back() = decomposition.dbegin;
	decomposition.dmap = new serializededata<esint, DIndex>(tarray<esint>(distribution, DOFDistribution), tarray<DIndex>(datadistribution, DOFData));
	esint index = 0;
	for (auto dmap = decomposition.dmap->cbegin(); dmap != decomposition.dmap->cend(); ++dmap, ++index) {
		if (dmap->size() > 1) {
			decomposition.sharedDOFs.push_back(index);
		}
	}

	decomposition.begin = dofs * info::mesh->nodes->uniqInfo.offset;
	decomposition.end = dofs * (info::mesh->nodes->uniqInfo.offset + info::mesh->nodes->uniqInfo.size);
	decomposition.totalSize = dofs * info::mesh->nodes->uniqInfo.totalSize;

	decomposition.neighbors = info::mesh->neighbors;
	decomposition.neighDOF.resize(decomposition.neighbors.size() + 1, decomposition.begin); // the last is my offset
	decomposition.halo.clear();
	decomposition.halo.reserve(dofs * info::mesh->nodes->uniqInfo.nhalo);
	for (esint n = 0; n < info::mesh->nodes->uniqInfo.nhalo; ++n) {
		for (int dof = 0; dof < dofs; ++dof) {
			decomposition.halo.push_back(dofs * info::mesh->nodes->uniqInfo.position[n] + dof);
		}
	}

	std::vector<esint> dBuffer = { decomposition.begin };
	if (!Communication::gatherUniformNeighbors(dBuffer, decomposition.neighDOF, decomposition.neighbors)) {
		eslog::internalFailure("cannot exchange matrix distribution info.\n");
	}
}

void buildPattern(UniformNodesFETIPattern *pattern, int dofs, DOFsDecomposition &decomposition, Matrix_Shape shape, int domain)
{
	auto sizeTriangle = [] (int size) {
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

	auto fillLower = [] (IJ* &target, const esint *begin, const esint *end) {
		for (auto row = begin, colbegin = begin; row != end; ++row, ++colbegin) {
			for (auto col = colbegin; col != end; ++col, ++target) {
				if (*row <= *col) {
					target->row = *col;
					target->column = *row;
				} else {
					target->row = *row;
					target->column = *col;
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

	switch (shape) {
	case Matrix_Shape::FULL: size = sizeFull; fill = fillFull; break;
	case Matrix_Shape::UPPER: size = sizeTriangle; fill = fillUpper; break;
	case Matrix_Shape::LOWER: size = sizeTriangle; fill = fillLower; break;
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
				auto DOFs = (decomposition.dmap->begin() + (*n * dofs + dof))->begin();
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

	pattern->bregion[domain].resize(info::mesh->boundaryRegions.size());
	std::vector<esint> belement(dofs * 8);
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			Ksize = 0; RHSsize = 0;
			for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
				auto element = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
				for (esint e = info::mesh->boundaryRegions[r]->eintervals[i].begin; e < info::mesh->boundaryRegions[r]->eintervals[i].end; ++e, ++element) {
					RHSsize += element->size() * dofs;
					Ksize += size(element->size() * dofs);
				}
			}
			pattern->bregion[domain][r].f.reserve(RHSsize);
			pattern->bregion[domain][r].K.reserve(Ksize);

			for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
				auto element = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
				for (esint e = info::mesh->boundaryRegions[r]->eintervals[i].begin; e < info::mesh->boundaryRegions[r]->eintervals[i].end; ++e, ++element) {
					belement.clear();
					for (int dof = 0; dof < dofs; ++dof) {
						for (auto n = element->begin(); n != element->end(); ++n) {
							auto DOFs = (decomposition.dmap->begin() + (*n * dofs + dof))->begin();
							while (DOFs->domain != domain + info::mesh->domains->offset) {
								++DOFs;
							}
							belement.push_back(DOFs->index);
						}
					}
					for (size_t i = 0; i < belement.size(); ++i) {
						pattern->bregion[domain][r].f.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), belement[i]) - RHSPattern.begin());
						for (size_t j = 0; j < belement.size(); ++j) {
							pattern->bregion[domain][r].K.push_back(std::lower_bound(KPattern.begin(), KPattern.end(), IJ{belement[i], belement[j]}) - KPattern.begin());
						}
					}
				}
			}
		} else {
			// TODO: how to set node regions?
//			pattern->bregion[r].f.reserve(info::mesh->boundaryRegions[r]->nodes->datatarray().size());
//			for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().cbegin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); ++n) {
//				belement.clear();
//				for (int dof = 0; dof < dofs; ++dof) {
//					belement.push_back(info::mesh->nodes->uniqInfo.position[*n] * dofs + dof);
//				}
//				for (size_t i = 0; i < belement.size(); ++i) {
//					pattern->bregion[r].f.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), belement[i]) - RHSPattern.begin());
//				}
//			}
		}
	}
}

static void dirichlet(UniformNodesFETIPattern *pattern, std::map<std::string, ECFExpression> &settings, int dofs)
{
	size_t size = 0;
	std::vector<esint> indices;
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (settings.find(region->name) != settings.end()) {
			size += region->nodes->datatarray().size();
			pattern->dirichletInfo.dirichlet = true;
		}
	}

	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (pattern->dirichletInfo.dirichlet) {
			for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
				for (int d = 0; d < dofs; ++d) {
					indices.push_back(*n * dofs + d);
				}
			}
		}
	}
	pattern->dirichletInfo.size = dofs * info::mesh->nodes->uniqInfo.offset;
	pattern->dirichletInfo.f = indices; // use the first region to store indices permutation;
	utils::sortAndRemoveDuplicates(indices);
	pattern->dirichletInfo.indices = indices;
	for (size_t i = 0; i < pattern->dirichletInfo.f.size(); ++i) {
		pattern->dirichletInfo.f[i] = std::lower_bound(indices.begin(), indices.end(), pattern->dirichletInfo.f[i]) - indices.begin();
	}
}

void UniformNodesFETIPattern::set(std::map<std::string, ECFExpression> &settings, int dofs, DOFsDecomposition &decomposition, Matrix_Shape shape)
{
	fillDecomposition(this, dofs, decomposition);
	for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
		buildPattern(this, dofs, decomposition, shape, domain);
	}
	dirichlet(this, settings, dofs);
}

void UniformNodesFETIPattern::fillCSR(esint *rows, esint *cols, esint domain)
{
	rows[0] = Indexing::CSR;
	cols[0] = elements[domain].column.front() + Indexing::CSR;
	size_t r = 1;
	for (size_t c = 1; c < elements[domain].column.size(); ++c) {
		cols[c] = elements[domain].column[c] + Indexing::CSR;
		if (elements[domain].row[c] != elements[domain].row[c - 1]) {
			rows[r++] = c + Indexing::CSR;
		}
	}
	rows[r] = elements[domain].column.size() + Indexing::CSR;
}
