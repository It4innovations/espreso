
#ifndef SRC_ANALYSIS_PATTERN_MATRIX_UNIFORM_DIRECT_H_
#define SRC_ANALYSIS_PATTERN_MATRIX_UNIFORM_DIRECT_H_

#include "apply.h"
#include "decomposition.direct.h"
#include "synchronization.h"
#include "analysis/math/matrix_distributed.h"
#include "math/primitives/matrix_info.h"
#include "esinfo/ecfinfo.h"

#include <vector>

namespace espreso {

struct MatrixUniformDirect {

    template <typename T>
    struct Sync: public Synchronization<Matrix_Distributed<T> > {
        std::vector<std::vector<T> > sBuffer, rBuffer;
        std::vector<std::vector<esint> > rOffset;
        std::vector<esint> nOffset;
        std::vector<int> neighbors;

        void init(MatrixUniformDirect &m);
        void gatherFromUpper(Matrix_Distributed<T> &m);
        void scatterToUpper(Matrix_Distributed<T> &m);
    };

    template <typename T>
    struct Apply: public ApplyMatrix<Matrix_Distributed<T>, T> {
        Matrix_CSR<T, esint> localM;
        Vector_Dense<T, esint> localV;
        std::vector<std::vector<T> > sBuffer, rBuffer;
        std::vector<std::vector<esint> > rOffset, sOffset;
        std::vector<int> neighbors;
        std::vector<esint> nDOF;
        esint offset;

        SpBLAS<Matrix_CSR, T, esint> spblas;

        void init(MatrixUniformDirect &m);
        void apply(Matrix_Distributed<T> &m, Vector_Distributed<Vector_Dense, T> &y, const T &alpha, const T &beta, const Vector_Distributed<Vector_Dense, T> &x);
    };

    MatrixUniformDirect(int DOFs);
    MatrixUniformDirect(HeatTransferLoadStepConfiguration &configuration, int multiplicity);
    MatrixUniformDirect(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity);

    template <typename T>
    void set(Matrix_Base<T> *m, Sync<T> &sync, Apply<T> &apply)
    {
        auto *_m = dynamic_cast<Matrix_Distributed<T> *>(m);
        _m->decomposition = &decomposition;
        _m->sync = &sync;
        _m->applyMatrix = &apply;
        _m->cluster.type = _m->type = type;
        _m->cluster.shape = _m->shape = Matrix_Shape::FULL; // always full
        // we set square matrix in order to be able to call local operations (e.g., apply)
        _m->cluster.resize(pattern.nrows, pattern.ncols, pattern.column.size());
        for (size_t r = 0; r < pattern.row.size(); ++r) {
            _m->cluster.rows[r] = pattern.row[r] + Indexing::CSR;
        }
        for (size_t c = 0; c < pattern.column.size(); ++c) {
            _m->cluster.cols[c] = pattern.column[c] + Indexing::CSR;
        }
    }

    template <typename T>
    void map(Matrix_Base<T> *m)
    {
        auto *_m = dynamic_cast<Matrix_Distributed<T> *>(m);
        _m->mapping.elements.resize(elements.offset.size());
        for (size_t i = 0; i < elements.offset.size(); ++i) {
            _m->mapping.elements[i].data = _m->cluster.vals;
            _m->mapping.elements[i].position = elements.permutation.data() + elements.offset[i];
        }

        _m->mapping.boundary.resize(boundary.size());
        for (size_t r = 0; r < boundary.size(); ++r) {
            _m->mapping.boundary[r].resize(boundary[r].offset.size());
            for (size_t i = 0; i < boundary[r].offset.size(); ++i) {
                _m->mapping.boundary[r][i].data = _m->cluster.vals;
                _m->mapping.boundary[r][i].position = boundary[r].permutation.data() + boundary[r].offset[i];
            }
        }
    }

    DecompositionDirect decomposition;

protected:
    struct RegionInfo {
        std::vector<esint> permutation;
        std::vector<esint> offset;
    };

    struct {
        esint nrows = 0, ncols = 0;
        std::vector<esint> row, column; // row, column indices
    } pattern;

    int dofs;
    Matrix_Type type;
    Matrix_Shape shape;
    RegionInfo elements;
    std::vector<RegionInfo> boundary; // RegionInfo per domain per boundary region

private:
    void buildPattern(int dofs);
};

}

#endif /* SRC_ANALYSIS_PATTERN_MATRIX_UNIFORM_DIRECT_H_ */
