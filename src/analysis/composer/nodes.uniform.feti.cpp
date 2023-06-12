
#include "nodes.uniform.feti.h"
#include "ij.h"

#include "basis/containers/serializededata.h"
#include "basis/containers/allocators.h"
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

	// fill DMAP
	decomposition.dmap = new serializededata<esint, DIndex>(
			tarray<esint>(info::mesh->nodes->domains->boundarytarray()),
			tarray<DIndex>(info::mesh->nodes->domains->datatarray().distribution(), 1));

	// order: inner, lambdas, dirichlet
	{ // go through the rest dofs
		esint index = 0;
		auto fix = decomposition.fixedDOFs.begin();
		auto dmap = decomposition.dmap->begin();
		for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains, ++dmap, ++index) {
			if (fix == decomposition.fixedDOFs.end() || *fix != index) {
				if (domains->size() == 1) {
					dmap->begin()->domain = *domains->begin();
					dmap->begin()->index = pattern->elements[dmap->begin()->domain - decomposition.dbegin].size++;
				}
			}
			if (fix != decomposition.fixedDOFs.end() && *fix == index) ++fix;
		}
	}
	{ // go through shared nodes
		esint index = 0;
		auto fix = decomposition.fixedDOFs.begin();
		auto dmap = decomposition.dmap->begin();
		for (auto domains = info::mesh->nodes->domains->begin(); domains != info::mesh->nodes->domains->end(); ++domains, ++dmap, ++index) {
			if (fix == decomposition.fixedDOFs.end() || *fix != index) {
				if (domains->size() > 1) {
					decomposition.sharedDOFs.push_back(index);
					auto di = dmap->begin();
					for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
						di->domain = *d;
						if (decomposition.ismy(*d)) {
							di->index = pattern->elements[*d - decomposition.dbegin].size++;
						}
					}
				}
			}
			if (fix != decomposition.fixedDOFs.end() && *fix == index) ++fix;
		}
	}
	// go through dirichlet
	for (auto i = decomposition.fixedDOFs.begin(); i != decomposition.fixedDOFs.end(); ++i) {
		auto domains = info::mesh->nodes->domains->begin() + *i;
		auto dmap = decomposition.dmap->begin() + *i;
		auto di = dmap->begin();
		for (auto d = domains->begin(); d != domains->end(); ++d, ++di) {
			di->domain = *d;
			if (decomposition.ismy(*d)) {
				di->index = pattern->elements[*d - decomposition.dbegin].size++;
			}
		}
	}
}

void buildPattern(UniformNodesFETIPattern *pattern, int dofs, DOFsDecomposition &decomposition, Matrix_Shape shape, int domain)
{
	auto ebegin = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain];
	auto eend = info::mesh->elements->nodes->cbegin() + info::mesh->domains->elements[domain + 1];

	size_t domainSize = 0, ni = 0;
	std::vector<esint> imap(info::mesh->nodes->size, -1);
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
			begin[imap[*n]] += enodes->size() - 1; // do not count diagonal
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
						++dataSize; indices[end[imap[*from]]++] = imap[*to];
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

	pattern->elements[domain].row.reserve(patternSize);
	pattern->elements[domain].column.reserve(patternSize);
	for (esint r = 0, ii = 0; r < domainSize; ++r) {
		for (esint c = begin[r]; c < end[r]; ++c) {
			pattern->elements[domain].row.push_back(r);
			pattern->elements[domain].column.push_back(indices[c]);
		}
	}

	pattern->elements[domain].K.reserve(dataSize);
	pattern->elements[domain].f.reserve(dataSize);
	std::vector<esint> local(20);
	for (auto enodes = ebegin; enodes != eend; ++enodes) {
		local.clear();
		for (auto from = enodes->begin(); from != enodes->end(); ++from) {
			if (shape == Matrix_Shape::FULL) {
				esint offset = begin[imap[*from]];
				for (auto to = enodes->begin(); to != enodes->end(); ++to) {
					esint coffset = 0;
					while (indices[offset + coffset] < imap[*to]) ++coffset;
					pattern->elements[domain].K.push_back(ioffset[imap[*from]] + coffset);
				}
			} else {
				for (auto to = from; to != enodes->end(); ++to) {
					esint min = std::min(imap[*from], imap[*to]), max = std::max(imap[*from], imap[*to]);
					esint offset = begin[min], coffset = 0;
					while (indices[offset + coffset] < max) ++coffset;
					pattern->elements[domain].K.push_back(ioffset[min] + coffset);
				}
			}
			pattern->elements[domain].f.push_back(imap[*from]);
		}
	}

	pattern->bregion[domain].resize(info::mesh->boundaryRegions.size());
	std::vector<esint> belement(dofs * 8);
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
				auto element = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
				for (esint e = info::mesh->boundaryRegions[r]->eintervals[i].begin; e < info::mesh->boundaryRegions[r]->eintervals[i].end; ++e, ++element) {
					for (auto from = element->begin(); from != element->end(); ++from) {
						pattern->bregion[domain][r].f.push_back(imap[*from]);
						if (shape == Matrix_Shape::FULL) {
							esint offset = begin[imap[*from]];
							for (auto to = element->begin(); to != element->end(); ++to) {
								esint coffset = 0;
								while (indices[offset + coffset] < imap[*to]) ++coffset;
								pattern->bregion[domain][r].K.push_back(ioffset[imap[*from]] + coffset);
							}
						} else {
							for (auto to = from; to != element->end(); ++to) {
								esint min = std::min(imap[*from], imap[*to]), max = std::max(imap[*from], imap[*to]);
								esint offset = begin[min], coffset = 0;
								while (indices[offset + coffset] < max) ++coffset;
								pattern->bregion[domain][r].K.push_back(ioffset[min] + coffset);
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

static void dirichlet(UniformNodesFETIPattern *pattern, std::map<std::string, ECFExpression> &settings, int dofs, DOFsDecomposition &decomposition)
{
	pattern->dirichletInfo.resize(info::mesh->boundaryRegions.size());
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (settings.find(region->name) != settings.end()) {
			pattern->dirichletInfo[r].dirichlet = true;
			for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
				for (int d = 0; d < dofs; ++d) {
					decomposition.fixedDOFs.push_back(*n * dofs + d);
				}
			}
		}
	}

	pattern->dirichletInfo[0].f = decomposition.fixedDOFs; // use the first region to store indices permutation;
	utils::sortAndRemoveDuplicates(decomposition.fixedDOFs);
	pattern->dirichletInfo[0].indices = decomposition.fixedDOFs;
	for (size_t i = 0; i < pattern->dirichletInfo[0].f.size(); ++i) {
		pattern->dirichletInfo[0].f[i] = std::lower_bound(decomposition.fixedDOFs.begin(), decomposition.fixedDOFs.end(), pattern->dirichletInfo[0].f[i]) - decomposition.fixedDOFs.begin();
	}
}

void UniformNodesFETIPattern::set(std::map<std::string, ECFExpression> &settings, int dofs, DOFsDecomposition &decomposition, Matrix_Shape shape)
{
	this->dofs = dofs;
	dirichlet(this, settings, dofs, decomposition);
	fillDecomposition(this, dofs, decomposition);
	#pragma omp parallel for
	for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
		buildPattern(this, dofs, decomposition, shape, domain);
	}
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
