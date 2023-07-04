
#ifndef SRC_ANALYSIS_BUILDER_UNIFORMBUILDER_DIRECT_H_
#define SRC_ANALYSIS_BUILDER_UNIFORMBUILDER_DIRECT_H_

#include "builder.h"
#include "basis/containers/serializededata.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "math/utils/distributed/distribution.h"
#include "math/utils/distributed/synchronization.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

struct UniformBuilderDirectPattern {

	struct RegionInfo {
		esint nrows = 0, ncols = 0, dirichlet = 0;
		std::vector<esint> row, column; // row, column indices
		std::vector<esint> A, b, indices; // local permutations
	};

	UniformBuilderDirectPattern(std::map<std::string, ECFExpression> &dirichlet, int dofs, Matrix_Shape shape);
	UniformBuilderDirectPattern(std::map<std::string, ECFExpressionOptionalVector> &dirichlet, int dofs, Matrix_Shape shape);
	~UniformBuilderDirectPattern();

protected:
	DOFsDistribution distribution;
	RegionInfo elements;
	std::vector<RegionInfo> bregion; // RegionInfo per domain per boundary region

private:
	void buildPattern(int dofs);
};

template <typename T>
struct UniformBuilderDirect: UniformBuilderDirectPattern, SparseMatrixBuilder<T> {

	UniformBuilderDirect(std::map<std::string, ECFExpression> &dirichlet, int dofs, Matrix_Shape shape)
	: UniformBuilderDirectPattern(dirichlet, dofs, shape), SparseMatrixBuilder<T>(dofs), syncM{}, syncV{}
	{

	}

	UniformBuilderDirect(std::map<std::string, ECFExpressionOptionalVector> &dirichlet, int dofs, Matrix_Shape shape)
	: UniformBuilderDirectPattern(dirichlet, dofs, shape), SparseMatrixBuilder<T>(dofs), syncM{}, syncV{}
	{

	}

	~UniformBuilderDirect()
	{
		if (syncV) { delete syncV; }
		if (syncM) { delete syncM; }
	}

	void fillDirichlet(Vector_Base<T> *v)
	{
		auto *_v = dynamic_cast<Vector_Distributed<Vector_Sparse, T> *>(v);
		_v->distribution = &distribution;
		_v->cluster.resize(elements.nrows, bregion[0].indices.size());
		for (size_t i = 0; i < bregion[0].indices.size(); ++i) {
			_v->cluster.indices[i] = bregion[0].indices[i];
		}
	}

	void fillVector(Vector_Base<T> *v)
	{
		auto *_v = dynamic_cast<Vector_Distributed<Vector_Dense, T> *>(v);
		_v->distribution = &distribution;
		_v->cluster.resize(elements.nrows);
		if (syncV == nullptr) {
			syncV = new Data_Synchronization<Vector_Dense, T>();
			syncV->init(*_v);
		}
		_v->synchronization = syncV;
	}

	void fillMatrix(Matrix_Base<T> *m, Matrix_Type type, Matrix_Shape shape)
	{
		auto *_m = dynamic_cast<Matrix_Distributed<Matrix_CSR, T> *>(m);
		_m->distribution = &distribution;
		_m->cluster.type = _m->type = type;
		_m->cluster.shape = _m->shape = Matrix_Shape::FULL; // always full
		// we set square matrix in order to be able to call local operations (e.g., apply)
		_m->cluster.resize(elements.nrows, elements.ncols, elements.row.size());
		_m->cluster.rows[0] = Indexing::CSR;
		_m->cluster.cols[0] = elements.column.front() + Indexing::CSR;
		size_t r = 1;
		for (size_t c = 1; c < elements.column.size(); ++c) {
			_m->cluster.cols[c] = elements.column[c] + Indexing::CSR;
			if (elements.row[c] != elements.row[c - 1]) {
				_m->cluster.rows[r++] = c + Indexing::CSR;
			}
		}
		_m->cluster.rows[r] = elements.column.size() + Indexing::CSR;

		if (syncM == nullptr) {
			syncM = new Data_Synchronization<Matrix_CSR, T>();
			syncM->init(*_m);
		}
		_m->synchronization = syncM;
	}

	void fillDirichletMap(Vector_Base<T> *v)
	{
		auto *_v = dynamic_cast<Vector_Distributed<Vector_Sparse, T> *>(v);
		_v->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 1, offset = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
			if (bregion[r].dirichlet) {
				int nsize = 0;
				for (int d = 0; d < this->dofs; ++d) {
					if (bregion[r].dirichlet & (1 << d)) {
						++nsize;
					}
				}
				_v->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->nodes->threads());
				for (size_t t = 0; t < _v->mapping.boundary[r].size(); ++t) {
					_v->mapping.boundary[r][t].filter = bregion[r].dirichlet;
					_v->mapping.boundary[r][t].data = _v->cluster.vals;
					_v->mapping.boundary[r][t].position = bregion[0].b.data() + offset;
					offset += region->nodes->datatarray().size(t) * nsize;
				}
			}
		}

	}

	void fillVectorMap(Vector_Base<T> *v)
	{
		auto *_v = dynamic_cast<Vector_Distributed<Vector_Dense, T> *>(v);
		_v->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (size_t i = 0, offset = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			_v->mapping.elements[i].data = _v->cluster.vals;
			_v->mapping.elements[i].position = elements.b.data() + offset;
			esint esize = this->dofs * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
			offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * esize;
		}

		_v->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				_v->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
				for (size_t i = 0, offset = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					_v->mapping.boundary[r][i].data = _v->cluster.vals;
					_v->mapping.boundary[r][i].position = bregion[r].b.data() + offset;
					esint esize = this->dofs * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
					offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * esize;
				}
			}
		}
	}

	void fillMatrixMap(Matrix_Base<T> *m)
	{
		auto *_m = dynamic_cast<Matrix_Distributed<Matrix_CSR, T> *>(m);
		_m->mapping.elements.resize(info::mesh->elements->eintervals.size());
		for (size_t i = 0, offset = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			_m->mapping.elements[i].data = _m->cluster.vals;
			_m->mapping.elements[i].position = elements.A.data() + offset;
			esint esize = this->dofs * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
			offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * esize * esize;
		}

		_m->mapping.boundary.resize(info::mesh->boundaryRegions.size());
		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			if (info::mesh->boundaryRegions[r]->dimension) {
				_m->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
				for (size_t i = 0, offset = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					_m->mapping.boundary[r][i].data = _m->cluster.vals;
					_m->mapping.boundary[r][i].position = bregion[r].A.data() + offset;
					esint esize = this->dofs * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
					offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * esize * esize;
				}
			}
		}
	}

protected:
	Data_Synchronization<Matrix_CSR, T> *syncM;
	Data_Synchronization<Vector_Dense, T> *syncV;
};

}

#endif /* SRC_ANALYSIS_BUILDER_UNIFORMBUILDER_DIRECT_H_ */
