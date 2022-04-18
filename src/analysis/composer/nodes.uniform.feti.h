
#ifndef SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_
#define SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_

#include "elementmapping.h"

#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "config/holders/expression.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "math/primitives/vector_dense.h"
#include "math/primitives/matrix_info.h"
#include "math/primitives/matrix_csr.h"
#include "math/generalization/vector_distributed.h"
#include "math/generalization/vector_feti.h"
#include "math/generalization/matrix_feti.h"
#include "math/utils/decomposed/decomposition.h"

namespace espreso {

struct UniformNodesFETIPattern {

	struct RegionInfo {
		esint size = 0, dirichlet = 0;
		std::vector<esint> row, column; // row, column indices
		std::vector<esint> K, f, indices; // local permutations
	};

	UniformNodesFETIPattern();
	~UniformNodesFETIPattern();

	void set(std::map<std::string, ECFExpression> &settings, int dofs, DOFsDecomposition &decomposition, Matrix_Shape shape);

	template<typename T>
	void fill(Vector_Distributed<Vector_Sparse, T> &v)
	{
		v.cluster.resize(dirichletInfo.size, dirichletInfo.indices.size());
		for (size_t i = 0; i < dirichletInfo.indices.size(); ++i) {
			v.cluster.indices[i] = dirichletInfo.indices[i];
		}
	}

	template<typename T>
	void fill(Vector_FETI<Vector_Dense, T> &v)
	{
		v.domains.resize(elements.size());
		#pragma omp parallel for
		for (size_t d = 0; d < v.domains.size(); ++d) {
			v.domains[d].resize(elements[d].size);
		}
	}

	template<typename T>
	void fill(Matrix_FETI<Matrix_CSR, T> &m)
	{
		m.domains.resize(elements.size());
		#pragma omp parallel for
		for (size_t d = 0; d < m.domains.size(); ++d) {
			m.domains[d].type = m.type;
			m.domains[d].shape = m.shape;
			m.domains[d].resize(elements[d].size, elements[d].size, elements[d].row.size());
			fillCSR(m.domains[d].rows, m.domains[d].cols, d);
		}
	}

	template<typename T>
	void setMap(Matrix_FETI<Matrix_CSR, T> *m) const
	{
		m->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
			for (esint i = info::mesh->elements->eintervalsDistribution[domain], offset = 0; i < info::mesh->elements->eintervalsDistribution[domain + 1]; ++i) {
				m->mapping.elements[i].data = m->domains[domain].vals;
				m->mapping.elements[i].position = elements[domain].K.data() + offset;
				offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
			}
		}
		m->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				m->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
				for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
					for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain], offset = 0; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
						m->mapping.boundary[r][i].data = m->domains[domain].vals;
						m->mapping.boundary[r][i].position = bregion[domain][r].K.data() + offset;
						offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
					}
				}
			}
		}
	}

	template<typename T>
	void setMap(Vector_FETI<Vector_Dense, T> *v) const
	{
		v->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
			for (esint i = info::mesh->elements->eintervalsDistribution[domain], offset = 0; i < info::mesh->elements->eintervalsDistribution[domain + 1]; ++i) {
				v->mapping.elements[i].data = v->domains[domain].vals;
				v->mapping.elements[i].position = elements[domain].f.data() + offset;
				offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
			}
		}

		v->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				v->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
				for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
					for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain], offset = 0; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
						v->mapping.boundary[r][i].data = v->domains[domain].vals;
						v->mapping.boundary[r][i].position = bregion[domain][r].f.data() + offset;
						offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
					}
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
			if (dirichletInfo.dirichlet) {
				v->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->nodes->threads());
				for (size_t t = 0; t < v->mapping.boundary[r].size(); ++t) {
					v->mapping.boundary[r][t].data = v->cluster.vals;
					v->mapping.boundary[r][t].position = dirichletInfo.f.data() + offset;
					offset += region->nodes->datatarray().size(t);
				}
			}
		}
	}

	void fillCSR(esint *rows, esint *cols, esint domain);

	int dofs;
	std::vector<RegionInfo> elements;
	std::vector<std::vector<RegionInfo> > bregion; // RegionInfo per domain per boundary region
	RegionInfo dirichletInfo;
};

}

#endif /* SRC_ANALYSIS_COMPOSER_NODES_UNIFORM_FETI_H_ */
