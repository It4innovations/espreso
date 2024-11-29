
#ifndef SRC_ANALYSIS_PATTERN_MATRIX_UNIFORM_FETI_H_
#define SRC_ANALYSIS_PATTERN_MATRIX_UNIFORM_FETI_H_

#include "apply.h"
#include "decomposition.feti.h"
#include "synchronization.h"
#include "analysis/math/matrix_feti.h"
#include "math/primitives/matrix_info.h"
#include "esinfo/ecfinfo.h"

#include <vector>

namespace espreso {

struct MatrixUniformFETI {

    template <typename T>
    struct Sync: public Synchronization<Matrix_FETI<T> > {
        void init(MatrixUniformFETI &m) {}
        void gatherFromUpper(Matrix_FETI<T> &m) {}
        void scatterToUpper(Matrix_FETI<T> &m) {}
    }; // no need for synchronization

    template <typename T>
    struct Apply: public ApplyMatrix<Matrix_FETI<T>, T> {
        std::vector<SpBLAS<Matrix_CSR, T, int> > spblas;
        Vector_FETI<Vector_Dense, T> in;
        struct {
            Vector_FETI<Vector_Dense, T> feti;
            Vector_Dense<T, esint> direct;
        } out;

        void init(MatrixUniformFETI &m);
        void apply(Matrix_FETI<T> &m, Vector_Distributed<Vector_Dense, T> &y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> &x);
    };

    MatrixUniformFETI(int DOFs);
    MatrixUniformFETI(HeatTransferLoadStepConfiguration &configuration, int multiplicity);
    MatrixUniformFETI(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity);

    template <typename T>
    void set(Matrix_Base<T> *m, Sync<T> &sync, Apply<T> &apply)
    {
        auto *_m = dynamic_cast<Matrix_FETI<T> *>(m);
        _m->decomposition = &decomposition;
//        _m->sync = &sync;
        _m->applyMatrix = &apply;
       _m->domains.resize(pattern.size());
       #pragma omp parallel for
       for (size_t d = 0; d < _m->domains.size(); ++d) {
           _m->domains[d].type = _m->type = type;
           _m->domains[d].shape = _m->shape = shape;
           if (BEM) {
               _m->domains[d].resize(pattern[d].surface, pattern[d].surface, pattern[d].surface * (pattern[d].surface - 1) / 2 + pattern[d].surface);
               _m->domains[d].rows[0] = Indexing::CSR;
               for (esint r = 0, cc = 0; r < pattern[d].surface; ++r) {
                   _m->domains[d].rows[r + 1] = _m->domains[d].rows[r] + pattern[d].surface - r;
                   for (esint c = r; c < pattern[d].surface; ++c) {
                       _m->domains[d].cols[cc++] = c + Indexing::CSR;
                   }
               }
           } else {
               _m->domains[d].resize(pattern[d].size, pattern[d].size, pattern[d].row.size());
               _m->domains[d].rows[0] = Indexing::CSR;
               _m->domains[d].cols[0] = pattern[d].column.front() + Indexing::CSR;
               size_t r = 1;
               for (size_t c = 1; c < pattern[d].column.size(); ++c) {
                   _m->domains[d].cols[c] = pattern[d].column[c] + Indexing::CSR;
                   if (pattern[d].row[c] != pattern[d].row[c - 1]) {
                       _m->domains[d].rows[r++] = c + Indexing::CSR;
                   }
               }
               _m->domains[d].rows[r] = pattern[d].column.size() + Indexing::CSR;
           }
       }
    }

    template <typename T>
    void map(Matrix_Base<T> *m)
    {
        auto *_m = dynamic_cast<Matrix_FETI<T> *>(m);
        size_t total = 0;
        for (size_t d = 0; d < elements.size(); ++d) {
            total += elements[d].offset.size();
        }
        _m->mapping.elements.resize(total);
        for (size_t d = 0, i = 0; d < elements.size(); ++d) {
            for (size_t j = 0; j < elements[d].offset.size(); ++j, ++i) {
                _m->mapping.elements[i].data = _m->domains[d].vals;
                _m->mapping.elements[i].position = elements[d].permutation.data() + elements[d].offset[j];
            }
        }

        // TODO
//        _m->mapping.boundary.resize(boundary.front().size());
//        for (size_t r = 0; r < boundary.front().size(); ++r) {
//            if () {
//                _m->mapping.boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());
//                for (esint domain = 0; domain < info::mesh->domains->size; ++domain) {
//                    for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[domain], offset = 0; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[domain + 1]; ++i) {
//                        _m->mapping.boundary[r][i].data = _m->domains[domain].vals;
//                        _m->mapping.boundary[r][i].position = bregion[domain][r].K.data() + offset;
//                        esint esize = this->dofs * Mesh::edata[info::mesh->boundaryRegions[r]->eintervals[i].code].nodes;
//                        switch (_m->domains[domain].shape) {
//                        case Matrix_Shape::FULL: esize = esize * esize; break;
//                        case Matrix_Shape::UPPER: esize = esize * (esize - 1) / 2 + esize; break;
//                        case Matrix_Shape::LOWER: esize = esize * (esize - 1) / 2 + esize; break;
//                        }
//                        offset += (info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin) * esize;
//                    }
//                }
//            }
//        }
    }

    DecompositionFETI decomposition;

protected:
    struct RegionInfo {
        std::vector<esint> permutation;
        std::vector<esint> offset;
    };

    struct Pattern {
        esint size = 0, surface = 0;
        std::vector<esint> row, column; // row, column indices
    };

    int dofs;
    bool BEM;
    Matrix_Type type;
    Matrix_Shape shape;
    std::vector<esint> surface;
    std::vector<Pattern> pattern;
    std::vector<RegionInfo> elements;
    std::vector<std::vector<RegionInfo> > boundary; // RegionInfo per domain per boundary region

private:
    void fillDecomposition(int dofs);
    void buildPattern(int dofs, Matrix_Shape shape, int domain);
};


}

#endif /* SRC_ANALYSIS_PATTERN_MATRIX_UNIFORM_FETI_H_ */
