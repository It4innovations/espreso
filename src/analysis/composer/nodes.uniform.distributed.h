
#ifndef SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_DISTRIBUTED_H_
#define SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_DISTRIBUTED_H_

#include "elementmapping.h"

#include "basis/containers/serializededata.h"
#include "basis/utilities/utils.h"
#include "esinfo/meshinfo.h"
#include "config/holders/expression.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "math2/primitives/vector_dense.h"
#include "math2/primitives/matrix_info.h"
#include "math2/primitives/matrix_csr.h"
#include "math2/generalization/vector_distributed.h"
#include "math2/generalization/matrix_distributed.h"
#include "math2/utils/dofs_distribution.h"
#include "math2/utils/utils_distributed.h"

#include <vector>

namespace espreso {

struct UniformNodesDistributedPattern {

	struct RegionInfo {
		esint nrows, ncols, dirichlet;
		std::vector<esint> row, column; // row, column indices
		std::vector<esint> A, b, indices; // local permutations
	};

	UniformNodesDistributedPattern();
	~UniformNodesDistributedPattern();

	void set(std::map<std::string, ECFExpression> &settings, int dofs, DOFsDistribution &distribution);

	template<typename T>
	void fill(Vector_Distributed<Vector_Sparse, T> &v)
	{
		v.cluster.resize(elements.nrows, bregion[0].indices.size());
		for (size_t i = 0; i < bregion[0].indices.size(); ++i) {
			v.cluster.indices[i] = bregion[0].indices[i];
		}
	}

	template<typename T>
	void fill(Vector_Distributed<Vector_Dense, T> &v)
	{
		v.cluster.resize(elements.nrows);
	}

	template<typename T>
	void fill(Matrix_Distributed<Matrix_CSR, T> &m)
	{
		m.cluster.type = m.type;
		m.cluster.shape = m.shape;
		// we set square matrix in order to be able to call local operations (e.g., apply)
		m.cluster.resize(elements.nrows, elements.ncols, elements.row.size());
		fillCSR(m.cluster.rows, m.cluster.cols);
	}

	template<typename T>
	void setMap(Matrix_Distributed<Matrix_CSR, T> *m) const
	{
		m->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (size_t i = 0, offset = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			m->mapping.elements[i].data = m->cluster.vals;
			m->mapping.elements[i].position = elements.A.data() + offset;
			offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
		}

		m->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				m->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
				for (size_t i = 0, offset = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					m->mapping.boundary[r][i].data = m->cluster.vals;
					m->mapping.boundary[r][i].position = bregion[r].A.data() + offset;
					offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
				}
			}
		}
	}

	template<typename T>
	void setMap(Vector_Distributed<Vector_Dense, T> *v) const
	{
		v->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (size_t i = 0, offset = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			v->mapping.elements[i].data = v->cluster.vals;
			v->mapping.elements[i].position = elements.b.data() + offset;
			offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
		}

		v->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				v->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
				for (size_t i = 0, offset = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					v->mapping.boundary[r][i].data = v->cluster.vals;
					v->mapping.boundary[r][i].position = bregion[r].b.data() + offset;
					offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
				}
			}
		}
	}

	template<typename T>
	void setDirichletMap(Vector_Distributed<Vector_Sparse, T> *v) const
	{
		v->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 1, offset = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
			if (bregion[r].dirichlet) {
				v->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->nodes->threads());
				for (size_t t = 0; t < v->mapping.boundary[r].size(); ++t) {
					v->mapping.boundary[r].front().data = v->cluster.vals;
					v->mapping.boundary[r].front().position = bregion[0].b.data() + offset;
					offset += region->nodes->datatarray().size(t);
				}
			}
		}
	}

	void fillCSR(esint *rows, esint *cols);

	int dofs;
	RegionInfo elements;
	std::vector<RegionInfo> bregion; // RegionInfo per domain per boundary region
};

}
#endif /* SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_DISTRIBUTED_H_ */
