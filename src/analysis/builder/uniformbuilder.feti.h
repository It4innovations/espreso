
#ifndef SRC_ANALYSIS_BUILDER_UNIFORMBUILDER_FETI_H_
#define SRC_ANALYSIS_BUILDER_UNIFORMBUILDER_FETI_H_

#include "builder.h"
#include "feti.decomposition.h"
#include "direct.synchronization.h"
#include "analysis/math/matrix_feti.h"
#include "analysis/math/vector_feti.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

struct UniformBuilderFETIPattern {

    struct RegionInfo {
        esint size = 0, dirichlet = 0;
        std::vector<esint> row, column; // row, column indices
        std::vector<esint> K, f, indices; // local permutations
    };

protected:
    int dofs;
    Matrix_Type type;
    Matrix_Shape shape;
    FETIDecomposition decomposition;
    std::vector<RegionInfo> elements;
    std::vector<std::vector<RegionInfo> > bregion; // RegionInfo per domain per boundary region
    std::vector<RegionInfo> dirichletInfo;

    UniformBuilderFETIPattern(HeatTransferLoadStepConfiguration &configuration);
    UniformBuilderFETIPattern(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity);

private:
    void fillDecomposition(FETIConfiguration &feti, int dofs);
    void buildPattern(int dofs, Matrix_Shape shape, int domain);
};

template <typename T>
struct UniformBuilderFETI: UniformBuilderFETIPattern, SparseMatrixBuilder<T> {

    UniformBuilderFETI(HeatTransferLoadStepConfiguration &configuration): UniformBuilderFETIPattern(configuration), syncS{} {}
    UniformBuilderFETI(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity): UniformBuilderFETIPattern(configuration, multiplicity), syncS{} {}

    ~UniformBuilderFETI()
    {
        if (syncS) { delete syncS; }
    }

    void fillDirichlet(Vector_Base<T> *v)
    {
        auto *_v = dynamic_cast<Vector_Distributed<Vector_Sparse, T> *>(v);
        _v->decomposition = &decomposition;
        _v->cluster.resize(dirichletInfo[0].size, dirichletInfo[0].indices.size());
        for (size_t i = 0; i < dirichletInfo[0].indices.size(); ++i) {
            _v->cluster.indices[i] = dirichletInfo[0].indices[i];
        }
        if (syncS == nullptr) {
            syncS = new Vector_Sparse_Sync<T>();
        }
        _v->_sync = syncS;
    }

    void fillVector(Vector_Base<T> *v)
    {
        auto *_v = dynamic_cast<Vector_FETI<Vector_Dense, T> *>(v);
        _v->decomposition = &decomposition;
        _v->domains.resize(elements.size());
        #pragma omp parallel for
        for (size_t d = 0; d < _v->domains.size(); ++d) {
            _v->domains[d].resize(elements[d].size);
        }
    }

    void fillMatrix(Matrix_Base<T> *m)
    {
        auto *_m = dynamic_cast<Matrix_FETI<T> *>(m);
        _m->decomposition = &decomposition;
        _m->domains.resize(elements.size());
        #pragma omp parallel for
        for (size_t d = 0; d < _m->domains.size(); ++d) {
            _m->domains[d].type = _m->type = type;
            _m->domains[d].shape = _m->shape = shape;
            _m->domains[d].resize(elements[d].size, elements[d].size, elements[d].row.size());
            _m->domains[d].rows[0] = Indexing::CSR;
            _m->domains[d].cols[0] = elements[d].column.front() + Indexing::CSR;
            size_t r = 1;
            for (size_t c = 1; c < elements[d].column.size(); ++c) {
                _m->domains[d].cols[c] = elements[d].column[c] + Indexing::CSR;
                if (elements[d].row[c] != elements[d].row[c - 1]) {
                    _m->domains[d].rows[r++] = c + Indexing::CSR;
                }
            }
            _m->domains[d].rows[r] = elements[d].column.size() + Indexing::CSR;
        }
    }

    void fillDirichletMap(Vector_Base<T> *v)
    {
        auto *_v = dynamic_cast<Vector_Distributed<Vector_Sparse, T> *>(v);
        _v->mapping.boundary.resize(info::mesh->boundaryRegions.size());
        for (size_t r = 1, offset = 0; r < info::mesh->boundaryRegions.size(); ++r) {
            const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
            if (dirichletInfo[r].dirichlet) {
                int nsize = 0;
                for (int d = 0; d < this->dofs; ++d) {
                    if (dirichletInfo[r].dirichlet & (1 << d)) {
                        ++nsize;
                    }
                }
                _v->mapping.boundary[r].resize(region->nodes->threads());
                for (size_t t = 0; t < _v->mapping.boundary[r].size(); ++t) {
                    _v->mapping.boundary[r][t].filter = dirichletInfo[r].dirichlet;
                    _v->mapping.boundary[r][t].data = _v->cluster.vals;
                    _v->mapping.boundary[r][t].position = dirichletInfo[0].f.data() + offset;
                    offset += region->nodes->datatarray().size(t) * nsize;
                }
            }
        }
    }

    void fillVectorMap(Vector_Base<T> *v)
    {
        auto *_v = dynamic_cast<Vector_FETI<Vector_Dense, T> *>(v);
        _v->mapping.elements.resize(info::mesh->elements->eintervals.size());
        for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
            for (esint i = info::mesh->elements->eintervalsDistribution[domain], offset = 0; i < info::mesh->elements->eintervalsDistribution[domain + 1]; ++i) {
                _v->mapping.elements[i].data = _v->domains[domain].vals;
                _v->mapping.elements[i].position = elements[domain].f.data() + offset;
                offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
            }
        }

        _v->mapping.boundary.resize(info::mesh->boundaryRegions.size());
        for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
            if (info::mesh->boundaryRegions[r]->dimension) {
                _v->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
                for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
                    for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain], offset = 0; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
                        _v->mapping.boundary[r][i].data = _v->domains[domain].vals;
                        _v->mapping.boundary[r][i].position = bregion[domain][r].f.data() + offset;
                        offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
                    }
                }
            }
        }
    }

    void fillMatrixMap(Matrix_Base<T> *m)
    {
        auto *_m = dynamic_cast<Matrix_FETI<T> *>(m);
        _m->mapping.elements.resize(info::mesh->elements->eintervals.size());
        for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
            for (esint i = info::mesh->elements->eintervalsDistribution[domain], offset = 0; i < info::mesh->elements->eintervalsDistribution[domain + 1]; ++i) {
                _m->mapping.elements[i].data = _m->domains[domain].vals;
                _m->mapping.elements[i].position = elements[domain].K.data() + offset;
                esint esize = this->dofs * Mesh::edata[info::mesh->elements->eintervals[i].code].nodes;
                switch (_m->domains[domain].shape) {
                case Matrix_Shape::FULL: esize = esize * esize; break;
                case Matrix_Shape::UPPER: esize = esize * (esize - 1) / 2 + esize; break;
                case Matrix_Shape::LOWER: esize = esize * (esize - 1) / 2 + esize; break;
                }
                offset += (info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin) * esize;
            }
        }
        _m->mapping.boundary.resize(info::mesh->boundaryRegions.size());
        for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
            if (info::mesh->boundaryRegions[r]->dimension) {
                _m->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
                for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
                    for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain], offset = 0; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
                        _m->mapping.boundary[r][i].data = _m->domains[domain].vals;
                        _m->mapping.boundary[r][i].position = bregion[domain][r].K.data() + offset;
                        esint esize = this->dofs * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
                        switch (_m->domains[domain].shape) {
                        case Matrix_Shape::FULL: esize = esize * esize; break;
                        case Matrix_Shape::UPPER: esize = esize * (esize - 1) / 2 + esize; break;
                        case Matrix_Shape::LOWER: esize = esize * (esize - 1) / 2 + esize; break;
                        }
                        offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * esize;
                    }
                }
            }
        }
    }

protected:
    Vector_Sparse_Sync<T> *syncS;
};

}

#endif /* SRC_ANALYSIS_BUILDER_UNIFORMBUILDER_FETI_H_ */
