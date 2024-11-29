
#include "analysis/math/matrix_distributed.h"
#include "analysis/math/vector_distributed.h"
#include "config/ecf/linearsolver/suitesparse.h"

namespace espreso {

template<typename T>
struct SuiteSparseDataHolder;

template<typename T>
class SuiteSparseDirectSolver {
public:
    SuiteSparseDirectSolver(SuiteSparseConfiguration &configuration)
    : configuration(configuration), external(nullptr)
    {
        check();
    }

    ~SuiteSparseDirectSolver()
    {
        clear();
    }

    bool set(const Matrix_Distributed<T> &A);
    bool update(const Matrix_Distributed<T> &A);
    bool solve(const Vector_Distributed<Vector_Dense, T> &b, Vector_Distributed<Vector_Dense, T> &x);
    bool solve(const Vector_Distributed<Matrix_Dense, T> &B, Vector_Distributed<Matrix_Dense, T> &X);

protected:
    SuiteSparseConfiguration &configuration;
    SuiteSparseDataHolder<T> *external;

    bool call(int phase);
    void check();
    void clear();
};

}
