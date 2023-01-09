
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
	dirichlet(this, settings, dofs, decomposition);
	fillDecomposition(this, dofs, decomposition);
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
