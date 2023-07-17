
#include "uniformbuilder.feti.h"

#include "basis/containers/serializededata.h"
#include "basis/containers/allocators.h"
#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "mesh/store/domainstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"
#include "mesh/store/nodestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

UniformBuilderFETIPattern::UniformBuilderFETIPattern(std::map<std::string, ECFExpression> &dirichlet, int dofs, Matrix_Shape shape)
{
	dirichletInfo.resize(info::mesh->boundaryRegions.size());
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (dirichlet.find(region->name) != dirichlet.end()) {
			dirichletInfo[r].dirichlet = 1;
		}
	}

	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (dirichlet.find(region->name) != dirichlet.end()) {
			if (dirichletInfo[r].dirichlet) {
				for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
					for (int d = 0; d < dofs; ++d) {
						decomposition.fixedDOFs.push_back(*n * dofs + d);
					}
				}
			}
		}
	}

	dirichletInfo[0].f = decomposition.fixedDOFs; // use the first region to store indices permutation;
	utils::sortAndRemoveDuplicates(decomposition.fixedDOFs);
	dirichletInfo[0].indices = decomposition.fixedDOFs;
	for (size_t i = 0; i < dirichletInfo[0].f.size(); ++i) {
		dirichletInfo[0].f[i] = std::lower_bound(decomposition.fixedDOFs.begin(), decomposition.fixedDOFs.end(), dirichletInfo[0].f[i]) - decomposition.fixedDOFs.begin();
	}

	fillDecomposition(dofs);
	#pragma omp parallel for
	for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
		buildPattern(dofs, shape, domain);
	}
}

UniformBuilderFETIPattern::UniformBuilderFETIPattern(std::map<std::string, ECFExpressionOptionalVector> &dirichlet, int dofs, Matrix_Shape shape)
{
	dirichletInfo.resize(info::mesh->boundaryRegions.size());
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		auto expr = dirichlet.find(region->name);
		if (expr != dirichlet.end()) {
			for (int d = 0; d < dofs; ++d) {
				if (expr->second.data[d].isset) {
					dirichletInfo[r].dirichlet += 1 << d;
				}
			}
		}
	}

	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (dirichlet.find(region->name) != dirichlet.end()) {
			for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
				for (int d = 0; d < dofs; ++d) {
					if (dirichletInfo[r].dirichlet & (1 << d)) {
						decomposition.fixedDOFs.push_back(*n * dofs + d);
					}
				}
			}
		}
	}

	dirichletInfo[0].f = decomposition.fixedDOFs; // use the first region to store indices permutation;
	utils::sortAndRemoveDuplicates(decomposition.fixedDOFs);
	dirichletInfo[0].indices = decomposition.fixedDOFs;
	for (size_t i = 0; i < dirichletInfo[0].f.size(); ++i) {
		dirichletInfo[0].f[i] = std::lower_bound(decomposition.fixedDOFs.begin(), decomposition.fixedDOFs.end(), dirichletInfo[0].f[i]) - decomposition.fixedDOFs.begin();
	}

	fillDecomposition(dofs);
	#pragma omp parallel for
	for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
		buildPattern(dofs, shape, domain);
	}
}

UniformBuilderFETIPattern::~UniformBuilderFETIPattern()
{

}

void UniformBuilderFETIPattern::fillDecomposition(int dofs)
{
	elements.resize(info::mesh->domains->size);
	bregion.resize(info::mesh->domains->size);

	decomposition.dbegin = info::mesh->domains->offset;
	decomposition.dend = info::mesh->domains->offset + info::mesh->domains->size;
	decomposition.dtotal = info::mesh->domains->totalSize;
	decomposition.neighDomain.resize(info::mesh->neighbors.size() + 1);
	auto ddist = info::mesh->domains->gatherProcDistribution(); // remove this
	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		decomposition.neighDomain[n] = ddist[info::mesh->neighbors[n]];
	}
	decomposition.neighDomain.back() = decomposition.dbegin;

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

	std::vector<esint> distribution(dofs * info::mesh->nodes->size + 1);
	auto domainmap = info::mesh->nodes->domains->begin();
	esint sum = 0;
	for (esint i = 0; i < info::mesh->nodes->size; ++i, ++domainmap) {
		for (int dof = 0; dof < dofs; ++dof, sum += domainmap->size()) {
			distribution[i * dofs + dof] = sum;
		}
	}
	distribution.back() = sum;

	// fill DMAP
	decomposition.dmap = new serializededata<esint, DIndex>(
			tarray<esint>(info::env::threads, distribution),
			tarray<DIndex>(info::mesh->nodes->domains->datatarray().distribution(), dofs));

	// order: inner, lambdas, dirichlet
	{ // go through the rest dofs
		esint index = 0;
		auto fix = decomposition.fixedDOFs.begin();
		auto dmap = decomposition.dmap->begin();
		for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains, ++index) {
			for (int dof = 0; dof < dofs; ++dof, ++dmap) {
				if (fix == decomposition.fixedDOFs.end() || *fix != dofs * index + dof) {
					if (domains->size() == 1) {
						dmap->begin()->domain = *domains->begin();
						dmap->begin()->index = elements[dmap->begin()->domain - decomposition.dbegin].size++;
					}
				}
				if (fix != decomposition.fixedDOFs.end() && *fix == dofs * index + dof) ++fix;
			}
		}
	}
	{ // go through shared nodes
		esint index = 0;
		auto fix = decomposition.fixedDOFs.begin();
		auto dmap = decomposition.dmap->begin();
		for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains, ++index) {
			for (int dof = 0; dof < dofs; ++dof, ++dmap) {
				if (fix == decomposition.fixedDOFs.end() || *fix != dofs * index + dof) {
					if (domains->size() > 1) {
						decomposition.sharedDOFs.push_back(dofs * index + dof);
						auto di = dmap->begin();
						for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
							di->domain = *d;
							if (decomposition.ismy(*d)) {
								di->index = elements[*d - decomposition.dbegin].size++;
							}
						}
					}
				}
				if (fix != decomposition.fixedDOFs.end() && *fix == dofs * index + dof) ++fix;
			}
		}
	}
	// go through dirichlet
	for (auto i = decomposition.fixedDOFs.begin(); i != decomposition.fixedDOFs.end(); ++i) {
		auto domains = info::mesh->nodes->domains->begin() + *i / dofs;
		auto dmap = decomposition.dmap->begin() + *i;
		auto di = dmap->begin();
		for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
			di->domain = *d;
			if (decomposition.ismy(*d)) {
				di->index = elements[*d - decomposition.dbegin].size++;
			}
		}
	}
}

void UniformBuilderFETIPattern::buildPattern(int dofs, Matrix_Shape shape, int domain)
{
	auto ebegin = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain];
	auto eend = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain + 1];

	size_t domainSize = 0, ni = 0;
	std::vector<esint> imap(dofs * info::mesh->nodes->size, -1);
	for (auto dmap = decomposition.dmap->cbegin(); dmap != decomposition.dmap->cend(); ++dmap, ++ni) {
		for (auto di = dmap->begin(); di != dmap->end(); ++di) {
			if (di->domain == domain + decomposition.dbegin) {
				++domainSize; imap[ni] = di->index;
			}
		}
	}

	std::vector<esint> begin(domainSize + 1, 1); // add diagonal
	for (auto enodes = ebegin; enodes != eend; ++enodes) {
		for (auto n = enodes->begin(); n != enodes->end(); ++n) {
			for (int dof = 0; dof < dofs; ++dof) {
				begin[imap[*n * dofs + dof]] += dofs * enodes->size() - 1; // do not count diagonal
			}
		}
	}
	utils::sizesToOffsets(begin);

	std::vector<esint> end = begin;
	std::vector<esint, initless_allocator<esint> > indices(begin.back());
	for (size_t n = 0; n < domainSize; ++n) {
		indices[end[n]++] = n; // diagonal
	}

	size_t dataSize = domainSize;
	for (auto enodes = ebegin; enodes != eend; ++enodes) {
		for (auto from = enodes->begin(); from != enodes->end(); ++from) {
			for (auto to = enodes->begin(); to != enodes->end(); ++to) {
				if (*from != *to) {
					if (shape == Matrix_Shape::FULL || imap[*from] < imap[*to]) {
						for (int d1 = 0; d1 < dofs; ++d1) {
							for (int d2 = 0; d2 < dofs; ++d2) {
								++dataSize; indices[end[imap[*from * dofs + d1]]++] = imap[*to * dofs + d2];
							}
						}
					}
				}
			}
		}
	}

	size_t patternSize = 0;
	std::vector<esint> ioffset(domainSize);
	for (esint n = 0; n < domainSize; ++n) {
		std::sort(indices.begin() + begin[n], indices.begin() + end[n]);
		esint unique = begin[n];
		for (auto i = begin[n] + 1; i < end[n]; ++i) {
			if (indices[unique] != indices[i]) {
				indices[++unique] = indices[i];
			}
		}
		end[n] = unique + 1;
		ioffset[n] = patternSize;
		patternSize += end[n] - begin[n];
	}

	elements[domain].row.reserve(patternSize);
	elements[domain].column.reserve(patternSize);
	for (esint r = 0, ii = 0; r < domainSize; ++r) {
		for (esint c = begin[r]; c < end[r]; ++c) {
			elements[domain].row.push_back(r);
			elements[domain].column.push_back(indices[c]);
		}
	}

	elements[domain].K.reserve(dataSize);
	elements[domain].f.reserve(dataSize);
	std::vector<esint> local(20);
	for (auto enodes = ebegin; enodes != eend; ++enodes) {
		local.clear();
		for (auto from = enodes->begin(); from != enodes->end(); ++from) {
			for (int d1 = 0; d1 < dofs; ++d1) {
				if (shape == Matrix_Shape::FULL) {
					esint offset = begin[imap[*from * dofs + d1]];
					for (auto to = enodes->begin(); to != enodes->end(); ++to) {
						for (int d2 = 0; d2 < dofs; ++d2) {
							esint coffset = 0;
							while (indices[offset + coffset] < imap[*to * dofs + d2]) ++coffset;
							elements[domain].K.push_back(ioffset[imap[*from * dofs + d1]] + coffset);
						}
					}
				} else {
					for (auto to = from; to != enodes->end(); ++to) {
						for (int d2 = 0; d2 < dofs; ++d2) {
							esint min = std::min(imap[*from * dofs + d1], imap[*to * dofs + d2]), max = std::max(imap[*from * dofs + d1], imap[*to * dofs + d2]);
							esint offset = begin[min], coffset = 0;
							while (indices[offset + coffset] < max) ++coffset;
							elements[domain].K.push_back(ioffset[min] + coffset);
						}
					}
				}
				elements[domain].f.push_back(imap[*from * dofs + d1]);
			}
		}
	}

	bregion[domain].resize(info::mesh->boundaryRegions.size());
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
				auto element = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
				for (esint e = info::mesh->boundaryRegions[r]->eintervals[i].begin; e < info::mesh->boundaryRegions[r]->eintervals[i].end; ++e, ++element) {
					for (auto from = element->begin(); from != element->end(); ++from) {
						for (int d1 = 0; d1 < dofs; ++d1) {
							bregion[domain][r].f.push_back(imap[*from * dofs + d1]);
							if (shape == Matrix_Shape::FULL) {
								esint offset = begin[imap[*from * dofs + d1]];
								for (auto to = element->begin(); to != element->end(); ++to) {
									for (int d2 = 0; d2 < dofs; ++d2) {
										esint coffset = 0;
										while (indices[offset + coffset] < imap[*to * dofs + d2]) ++coffset;
										bregion[domain][r].K.push_back(ioffset[imap[*from * dofs + d1]] + coffset);
									}
								}
							} else {
								for (auto to = from; to != element->end(); ++to) {
									for (int d2 = 0; d2 < dofs; ++d2) {
										esint min = std::min(imap[*from * dofs + d1], imap[*to * dofs + d2]), max = std::max(imap[*from * dofs + d1], imap[*to * dofs + d2]);
										esint offset = begin[min], coffset = 0;
										while (indices[offset + coffset] < max) ++coffset;
										bregion[domain][r].K.push_back(ioffset[min] + coffset);
									}
								}
							}
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

