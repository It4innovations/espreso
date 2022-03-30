
#include "nodes.uniform.distributed.h"
#include "ij.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/nodestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

UniformNodesDistributedPattern::UniformNodesDistributedPattern()
: dofs(0)
{

}

UniformNodesDistributedPattern::~UniformNodesDistributedPattern()
{

}

void buildPattern(UniformNodesDistributedPattern *pattern, int dofs, DOFsDistribution &distribution)
{
	double start = eslog::time();
	eslog::info(" == LINEAR SYSTEM                                                               DISTRIBUTED == \n");
	eslog::info(" == DOFS PER NODE                                                                         %d == \n", dofs);

	pattern->dofs = dofs;
	auto size = [] (int size) {
		return size * size;
	};

	auto fill = [] (IJ* &target, const esint *begin, const esint *end) {
		for (auto row = begin; row != end; ++row) {
			for (auto col = begin; col != end; ++col, ++target) {
				target->row = *row;
				target->column = *col;
			}
		}
	};

	esint Asize = 0, RHSsize = 0;
	for (auto e = info::mesh->elements->nodes->cbegin(); e != info::mesh->elements->nodes->cend(); ++e) {
		RHSsize += e->size() * dofs;
		Asize += size(e->size() * dofs);
	}

	std::vector<IJ> APattern(Asize);
	std::vector<esint> RHSPattern(RHSsize);

	IJ *Aoffset = APattern.data();
	esint *RHSoffset = RHSPattern.data();

	for (auto e = info::mesh->elements->nodes->cbegin(); e != info::mesh->elements->nodes->cend(); ++e) {
		esint *_RHS = RHSoffset;
		for (int dof = 0; dof < dofs; ++dof) {
			for (auto n = e->begin(); n != e->end(); ++n, ++RHSoffset) {
				*RHSoffset = info::mesh->nodes->uniqInfo.position[*n] * dofs + dof;
			}
		}
		fill(Aoffset, _RHS, RHSoffset);
	}

	std::vector<IJ> AData = APattern;
	std::vector<esint> RHSData = RHSPattern;

	utils::sortAndRemoveDuplicates(APattern);
	utils::sortAndRemoveDuplicates(RHSPattern);

	// send halo rows to the holder process
	std::vector<std::vector<IJ> > sPatter(info::mesh->neighbors.size()), rPatter(info::mesh->neighbors.size());
	size_t max_halo = std::lower_bound(APattern.begin(), APattern.end(), dofs * info::mesh->nodes->uniqInfo.nhalo + 1, [] (const IJ &ij, esint i) { return ij.row < i; }) - APattern.begin();
	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		sPatter[n].reserve(max_halo);
	}

	auto ranks = info::mesh->nodes->ranks->cbegin();
	auto begin = APattern.begin(), end = begin;
	for (esint n = 0; n < info::mesh->nodes->uniqInfo.nhalo; ++n, ++ranks, begin = end) {
		while (begin->row + dofs > (++end)->row);
		auto neigh = info::mesh->neighbors.begin();
		for (auto r = ranks->begin(); r != ranks->end(); ++r) {
			if (*r != info::mpi::rank) {
				while (*neigh < *r) { ++neigh; }
				sPatter[neigh - info::mesh->neighbors.begin()].insert(sPatter[neigh - info::mesh->neighbors.begin()].end(), begin, end);
			}
		}
	}

	if (!Communication::receiveUpperUnknownSize(sPatter, rPatter, info::mesh->neighbors)) {
		eslog::internalFailure("exchange K pattern.\n");
	}

	for (size_t n = 0; n < info::mesh->neighbors.size(); ++n) {
		APattern.insert(APattern.end(), rPatter[n].begin(), rPatter[n].end());
	}
	utils::sortAndRemoveDuplicates(APattern);

	pattern->elements.nrows = dofs * (info::mesh->nodes->uniqInfo.nhalo + info::mesh->nodes->uniqInfo.size);
	pattern->elements.ncols = dofs * info::mesh->nodes->uniqInfo.totalSize;
	pattern->elements.row.reserve(APattern.size());
	pattern->elements.column.reserve(APattern.size());
	pattern->elements.A.reserve(AData.size());
	pattern->elements.b.reserve(RHSData.size());

	for (size_t i = 0; i < APattern.size(); ++i) {
		pattern->elements.row.push_back(APattern[i].row);
		pattern->elements.column.push_back(APattern[i].column);
	}

	for (size_t i = 0; i < AData.size(); ++i) {
		pattern->elements.A.push_back(std::lower_bound(APattern.begin(), APattern.end(), AData[i]) - APattern.begin());
	}
	for (size_t i = 0; i < RHSData.size(); ++i) {
		pattern->elements.b.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), RHSData[i]) - RHSPattern.begin());
	}

	std::vector<esint> belement(dofs * 8);
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {

			Asize = 0; RHSsize = 0;
			for (auto e = info::mesh->boundaryRegions[r]->elements->cbegin(); e != info::mesh->boundaryRegions[r]->elements->cend(); ++e) {
				RHSsize += e->size() * dofs;
				Asize += size(e->size() * dofs);
			}
			pattern->bregion[r].b.reserve(RHSsize);
			pattern->bregion[r].A.reserve(Asize);

			for (auto e = info::mesh->boundaryRegions[r]->elements->cbegin(); e != info::mesh->boundaryRegions[r]->elements->cend(); ++e) {
				belement.clear();
				for (int dof = 0; dof < dofs; ++dof) {
					for (auto n = e->begin(); n != e->end(); ++n) {
						belement.push_back(info::mesh->nodes->uniqInfo.position[*n] * dofs + dof);
					}
				}
				for (size_t i = 0; i < belement.size(); ++i) {
					pattern->bregion[r].b.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), belement[i]) - RHSPattern.begin());
					for (size_t j = 0; j < belement.size(); ++j) {
						pattern->bregion[r].A.push_back(std::lower_bound(APattern.begin(), APattern.end(), IJ{belement[i], belement[j]}) - APattern.begin());
					}
				}
			}
		} else {
			pattern->bregion[r].b.reserve(info::mesh->boundaryRegions[r]->nodes->datatarray().size());
			for (auto n = info::mesh->boundaryRegions[r]->nodes->datatarray().cbegin(); n != info::mesh->boundaryRegions[r]->nodes->datatarray().end(); ++n) {
				belement.clear();
				for (int dof = 0; dof < dofs; ++dof) {
					belement.push_back(info::mesh->nodes->uniqInfo.position[*n] * dofs + dof);
				}
				for (size_t i = 0; i < belement.size(); ++i) {
					pattern->bregion[r].b.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), belement[i]) - RHSPattern.begin());
				}
			}
		}
	}

	distribution.begin = dofs * info::mesh->nodes->uniqInfo.offset;
	distribution.end = dofs * (info::mesh->nodes->uniqInfo.offset + info::mesh->nodes->uniqInfo.size);
	distribution.totalSize = dofs * info::mesh->nodes->uniqInfo.totalSize;

	distribution.neighbors = info::mesh->neighbors;
	distribution.neighDOF.resize(distribution.neighbors.size() + 1, distribution.begin); // the last is my offset
	distribution.halo.clear();
	distribution.halo.reserve(dofs * info::mesh->nodes->uniqInfo.nhalo);
	for (esint n = 0; n < info::mesh->nodes->uniqInfo.nhalo; ++n) {
		for (int dof = 0; dof < dofs; ++dof) {
			distribution.halo.push_back(dofs * info::mesh->nodes->uniqInfo.position[n] + dof);
		}
	}

	std::vector<esint> dBuffer = { distribution.begin };
	if (!Communication::gatherUniformNeighbors(dBuffer, distribution.neighDOF, distribution.neighbors)) {
		eslog::internalFailure("cannot exchange matrix distribution info.\n");
	}

	size_t nonzeros[2] = { pattern->elements.row.size(), pattern->bregion[0].b.size() };
	Communication::allReduce(&nonzeros, NULL, 2, MPITools::getType<size_t>().mpitype, MPI_SUM);

	eslog::info(" == DIRICHLET SIZE                                                           %14lu == \n", nonzeros[1]);
	eslog::info(" == LINEAR SYSTEM SIZE                                                       %14d == \n", distribution.totalSize);
	eslog::info(" == NON-ZERO VALUES                                                          %14lu == \n", nonzeros[0]);
	eslog::info(" == NON-ZERO FILL-IN RATIO                                                         %7.4f\% == \n", 100.0 * nonzeros[0] / distribution.totalSize / distribution.totalSize);
	eslog::info(" == COMPOSITION RUNTIME                                                          %8.3f s == \n", eslog::time() - start);
	eslog::info(" ============================================================================================= \n");
}

static void dirichlet(UniformNodesDistributedPattern *pattern, std::map<std::string, ECFExpression> &settings, int dofs)
{
	size_t size = 0;
	std::vector<esint> indices;
	pattern->bregion.resize(info::mesh->boundaryRegions.size());
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (settings.find(region->name) != settings.end()) {
			size += region->nodes->datatarray().size();
			pattern->bregion[r].dirichlet = true;
		}
	}

	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (pattern->bregion[r].dirichlet) {
			for (auto n = region->nodes->datatarray().cbegin(); n != region->nodes->datatarray().cend(); ++n) {
				for (int d = 0; d < dofs; ++d) {
					indices.push_back(*n * dofs + d);
				}
			}
		}
	}
	pattern->bregion[0].b = indices; // use the first region to store indices permutation;
	utils::sortAndRemoveDuplicates(indices);
	pattern->bregion[0].indices = indices;
	for (size_t i = 0; i < pattern->bregion[0].b.size(); ++i) {
		pattern->bregion[0].b[i] = std::lower_bound(indices.begin(), indices.end(), pattern->bregion[0].b[i]) - indices.begin();
	}
}

void UniformNodesDistributedPattern::set(std::map<std::string, ECFExpression> &settings, int dofs, DOFsDistribution &distribution)
{
	dirichlet(this, settings, dofs);
	buildPattern(this, dofs, distribution);
}

void UniformNodesDistributedPattern::fillCSR(esint *rows, esint *cols)
{
	rows[0] = Indexing::CSR;
	cols[0] = elements.column.front() + Indexing::CSR;
	size_t r = 1;
	for (size_t c = 1; c < elements.column.size(); ++c) {
		cols[c] = elements.column[c] + Indexing::CSR;
		if (elements.row[c] != elements.row[c - 1]) {
			rows[r++] = c + Indexing::CSR;
		}
	}
	rows[r] = elements.column.size() + Indexing::CSR;
}
