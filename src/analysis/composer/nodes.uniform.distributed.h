
#ifndef SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_DISTRIBUTED_H_
#define SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_DISTRIBUTED_H_

#include "elementmapping.h"

#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

#include "math2/primitives/vector_dense.h"
#include "math2/primitives/matrix_info.h"
#include "math2/primitives/matrix_csr.h"
#include "math2/generalization/vector_distributed.h"
#include "math2/generalization/matrix_distributed.h"
#include "math2/utils/dofs_distribution.h"

#include <vector>

namespace espreso {

struct UniformNodesDistributedPattern {

	struct RegionInfo {
		esint size;
		std::vector<esint> row, column; // row, column indices
		std::vector<esint> K, f; // permutations
	};

	UniformNodesDistributedPattern();
	~UniformNodesDistributedPattern();

	void set(int dofs);

	template<typename T>
	void fill(Vector_Distributed<Vector_Dense, T> &v)
	{
		v.cluster.resize(elements.size);
		fillDistribution(v.distribution);
	}

	template<typename T>
	void fill(Matrix_Distributed<Matrix_CSR, T> &m)
	{
		m.cluster.type = m.type;
		m.cluster.shape = m.shape;
		m.cluster.resize(elements.size, elements.size, elements.K.size());
		fillCSR(m.cluster.rows, m.cluster.cols);
		fillDistribution(m.distribution);
	}

	template<typename T>
	void setMap(Matrix_Distributed<Matrix_CSR, T> *m) const
	{
		m->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (size_t i = 0, offset = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			m->mapping.elements[i].data = m->cluster.vals;
			m->mapping.elements[i].position = elements.K.data() + offset;
			offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
		}
	}

	template<typename T>
	void setMap(Vector_Distributed<Vector_Dense, T> *v) const
	{
		v->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (size_t i = 0, offset = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			v->mapping.elements[i].data = v->cluster.vals + offset;
			v->mapping.elements[i].position = elements.f.data() + offset;
			offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
		}
	}

	void fillCSR(esint *rows, esint *cols);
	void fillDistribution(DOFsDistribution &distribution);

	int dofs;
	RegionInfo elements;
	std::vector<RegionInfo> bregion; // RegionInfo per domain per boundary region
};

}
#endif /* SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_DISTRIBUTED_H_ */
