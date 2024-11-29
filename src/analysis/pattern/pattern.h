
#ifndef SRC_ANALYSIS_PATTERN_PATTERN_H_
#define SRC_ANALYSIS_PATTERN_PATTERN_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_distributed.h"
#include "analysis/linearsystem/linearsystem.h"

#include "matrix.uniform.direct.h"
#include "matrix.uniform.feti.h"
#include "vector.uniform.dense.h"
#include "vector.uniform.sparse.h"


namespace espreso {

template <typename T>
struct Pattern {
    virtual ~Pattern() {};
    virtual void set(Matrix_Base<T> *m)                                =0;
    virtual void set(Vector_Distributed<Vector_Dense, T> *v)           =0;
    virtual void set(Vector_Distributed<Matrix_Dense, T> *v, int nrhs) =0;
    virtual void set(Vector_Distributed<Vector_Sparse, T> *v)          =0;
    virtual void set(LinearSystemSolver<T> *solver)                    =0;

    virtual void map(Matrix_Base<T> *m)                       =0;
    virtual void map(Vector_Distributed<Vector_Dense, T> *v)  =0;
    virtual void map(Vector_Distributed<Matrix_Dense, T> *v)  =0;
    virtual void map(Vector_Distributed<Vector_Sparse, T> *v) =0;
};

template <typename Matrix, typename VectorDense, typename VectorSparse, typename T>
struct PatternInstance: public Pattern<T> {

    PatternInstance(int DOFs): matrix(DOFs), vector(DOFs), dirichlet(DOFs)
    {
        sync.matrix.init(matrix);
        sync.denseVector.init(matrix.decomposition);
        sync.denseMatrix.init(matrix.decomposition);
        apply.matrix.init(matrix);
    }

    template <typename Configuration>
    PatternInstance(Configuration &configuration, int multiplicity): matrix(configuration, multiplicity), vector(configuration, multiplicity), dirichlet(configuration, multiplicity)
    {
        sync.matrix.init(matrix);
        sync.denseVector.init(matrix.decomposition);
        apply.matrix.init(matrix);
    }

    void set(Matrix_Base<T> *m)
    {
        matrix.set(m, sync.matrix, apply.matrix);
    }

    void set(Vector_Distributed<Vector_Dense, T> *v)
    {
        vector.set(matrix.decomposition, v);
        v->sync = &sync.denseVector;
    }

    void set(Vector_Distributed<Matrix_Dense, T> *v, int nrhs)
    {
        vector.set(matrix.decomposition, v, nrhs);
        v->sync = &sync.denseMatrix;
    }

    void set(Vector_Distributed<Vector_Sparse, T> *v)
    {
        dirichlet.set(matrix.decomposition, v);
        v->sync = &sync.dirichlet;
    }

    void set(LinearSystemSolver<T> *solver)
    {
        set(solver->A);
        set(solver->b);
        set(solver->x);
        set(solver->dirichlet);
    }

    void map(Matrix_Base<T> *m)
    {
        matrix.map(m);
    }

    void map(Vector_Distributed<Vector_Dense, T> *v)
    {
        vector.map(v);
    }

    void map(Vector_Distributed<Matrix_Dense, T> *v)
    {
        vector.map(v);
    }

    void map(Vector_Distributed<Vector_Sparse, T> *v)
    {
        dirichlet.map(v);
    }

    Matrix matrix;
    VectorDense vector;
    VectorSparse dirichlet;

    struct {
        typename Matrix::template Sync<T> matrix;
        typename VectorDense::template SyncVector<T> denseVector;
        typename VectorDense::template SyncMatrix<T> denseMatrix;
        typename VectorSparse::template Sync<T> dirichlet;
    } sync;

    struct {
        typename Matrix::template Apply<T> matrix;
    } apply;
};

template <typename T> struct PatternUniformDirect: public PatternInstance<MatrixUniformDirect, VectorUniformDense, VectorUniformSparse, T> {
    PatternUniformDirect(int DOFs)                                                                 : PatternInstance<MatrixUniformDirect, VectorUniformDense, VectorUniformSparse, T>(DOFs) {}
    PatternUniformDirect(HeatTransferLoadStepConfiguration        &configuration, int multiplicity): PatternInstance<MatrixUniformDirect, VectorUniformDense, VectorUniformSparse, T>(configuration, multiplicity) {}
    PatternUniformDirect(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity): PatternInstance<MatrixUniformDirect, VectorUniformDense, VectorUniformSparse, T>(configuration, multiplicity) {}
};

template <typename T> struct PatternUniformFETI:   public PatternInstance<MatrixUniformFETI  , VectorUniformDense, VectorUniformSparse, T> {
    PatternUniformFETI(int DOFs)                                                                 : PatternInstance<MatrixUniformFETI, VectorUniformDense, VectorUniformSparse, T>(DOFs) {}
    PatternUniformFETI(HeatTransferLoadStepConfiguration        &configuration, int multiplicity): PatternInstance<MatrixUniformFETI, VectorUniformDense, VectorUniformSparse, T>(configuration, multiplicity) {}
    PatternUniformFETI(StructuralMechanicsLoadStepConfiguration &configuration, int multiplicity): PatternInstance<MatrixUniformFETI, VectorUniformDense, VectorUniformSparse, T>(configuration, multiplicity) {}
};

}

#endif /* SRC_ANALYSIS_PATTERN_PATTERN_H_ */
