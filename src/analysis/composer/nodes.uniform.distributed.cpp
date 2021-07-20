
#include "nodes.uniform.distributed.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"

#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/nodestore.h"

#include "wrappers/mpi/communication.h"

using namespace espreso;

struct IJ { esint row, column; };

inline bool operator==(const IJ &left, const IJ &right) { return left.row == right.row && left.column == right.column; }
inline bool operator!=(const IJ &left, const IJ &right) { return !(left == right); }
inline bool operator <(const IJ &left, const IJ &right) { return left.row == right.row ? left.column < right.column : left.row < right.row; }

UniformNodesDistributedPattern::UniformNodesDistributedPattern()
: dofs(0)
{

}

UniformNodesDistributedPattern::~UniformNodesDistributedPattern()
{

}

void fillPermutation(UniformNodesDistributedPattern *pattern, int dofs)
{
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

	esint Ksize = 0, RHSsize = 0;
	for (auto e = info::mesh->elements->nodes->cbegin(); e != info::mesh->elements->nodes->cend(); ++e) {
		RHSsize += e->size() * dofs;
		Ksize += size(e->size() * dofs);
	}

	std::vector<IJ> KPattern(Ksize);
	std::vector<esint> RHSPattern(RHSsize);

	IJ *Koffset = KPattern.data();
	esint *RHSoffset = RHSPattern.data();

	for (auto e = info::mesh->elements->nodes->cbegin(); e != info::mesh->elements->nodes->cend(); ++e) {
		esint *_RHS = RHSoffset;
		for (int dof = 0; dof < dofs; ++dof) {
			for (auto n = e->begin(); n != e->end(); ++n, ++RHSoffset) {
				*RHSoffset = info::mesh->nodes->uniqInfo.position[*n] * dofs + dof;
			}
		}
		fill(Koffset, _RHS, RHSoffset);
	}

	std::vector<IJ> KData = KPattern;
	std::vector<esint> RHSData = RHSPattern;

	utils::sortAndRemoveDuplicates(KPattern);
	utils::sortAndRemoveDuplicates(RHSPattern);

	pattern->elements.size = info::mesh->nodes->uniqInfo.nhalo + info::mesh->nodes->uniqInfo.size;
	pattern->elements.row.reserve(KPattern.size());
	pattern->elements.column.reserve(KPattern.size());
	pattern->elements.K.reserve(KData.size());
	pattern->elements.f.reserve(RHSData.size());

	for (size_t i = 0; i < KPattern.size(); ++i) {
		pattern->elements.row.push_back(KPattern[i].row);
		pattern->elements.column.push_back(KPattern[i].column);
	}

	for (size_t i = 0; i < KData.size(); ++i) {
		pattern->elements.K.push_back(std::lower_bound(KPattern.begin(), KPattern.end(), KData[i]) - KPattern.begin());
	}
	for (size_t i = 0; i < RHSData.size(); ++i) {
		pattern->elements.f.push_back(std::lower_bound(RHSPattern.begin(), RHSPattern.end(), RHSData[i]) - RHSPattern.begin());
	}
}

void UniformNodesDistributedPattern::set(int dofs)
{
	fillPermutation(this, dofs);
}

void UniformNodesDistributedPattern::fillCSR(esint *rows, esint *cols)
{
	rows[0] = cols[0] = _Matrix_CSR_Pattern::Indexing;
	size_t r = 1;
	for (size_t c = 1; c < elements.column.size(); ++c) {
		cols[c] = elements.column[c] + _Matrix_CSR_Pattern::Indexing;
		if (elements.row[c] != elements.row[c - 1]) {
			rows[r++] = c + _Matrix_CSR_Pattern::Indexing;
		}
	}
	rows[r] = elements.column.size() + _Matrix_CSR_Pattern::Indexing;
}

void UniformNodesDistributedPattern::fillDistribution(DOFsDistribution &distribution)
{
	distribution.begin = dofs * info::mesh->nodes->uniqInfo.offset;
	distribution.end = dofs * (info::mesh->nodes->uniqInfo.offset + info::mesh->nodes->uniqInfo.size);
	distribution.totalSize = dofs * info::mesh->nodes->uniqInfo.totalSize;

	distribution.neighbors = info::mesh->neighbors;
	distribution.neighDOF.resize(distribution.neighbors.size());
	distribution.halo.clear();
	distribution.halo.reserve(dofs * info::mesh->nodes->uniqInfo.nhalo);
	for (esint n = 0; n < info::mesh->nodes->uniqInfo.nhalo; ++n) {
		for (int dof = 0; dof < dofs; ++dof) {
			distribution.halo.push_back(dofs * info::mesh->nodes->uniqInfo.position[n] + dof);
		}
	}

	std::vector<esint> sBuffer = { distribution.begin };
	if (!Communication::gatherUniformNeighbors(sBuffer, distribution.neighDOF, distribution.neighbors)) {
		eslog::internalFailure("cannot exchange matrix distribution info.\n");
	}
}
