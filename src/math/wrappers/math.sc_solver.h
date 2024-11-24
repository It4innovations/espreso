
#ifndef SRC_MATH_WRAPPERS_MATH_SC_SOLVER_H_
#define SRC_MATH_WRAPPERS_MATH_SC_SOLVER_H_

#include "math/math.h"



namespace espreso {

template<typename T, typename I>
struct Schur_Complement_Solver_External_Representation;

template <typename T, typename I = int>
class SchurComplementSolver {

public:

    SchurComplementSolver();
    SchurComplementSolver(const SchurComplementSolver & other) = delete;
    SchurComplementSolver(SchurComplementSolver && other);
    SchurComplementSolver & operator=(const SchurComplementSolver & other) = delete;
    SchurComplementSolver & operator=(SchurComplementSolver && other);
    ~SchurComplementSolver();

    static const char * name();

    // stage = 0
    void commitMatrix(const Matrix_CSR<T,I> & A, I sc_size);
    void commitMatrix(const Matrix_CSR<T,I> & A11, const Matrix_CSR<T,I> & A12, const Matrix_CSR<T,I> & A21, const Matrix_CSR<T,I> & A22);
    // stage = 1
    void factorizeSymbolic();
    // stage = 2
    void updateMatrixValues();
    // stage = 3
    template<typename A>
    void factorizeNumericAndGetSc(Matrix_Dense<T,I,A> & sc, char uplo, T alpha = T{1});
    // stage = 4
    void solveA11(const Vector_Dense<T,I> & rhs, Vector_Dense<T,I> & sol);

    void say_hello();

private:
    std::unique_ptr<Schur_Complement_Solver_External_Representation<T,I>> ext;
};

}

#endif /* SRC_MATH_WRAPPERS_MATH_SC_SOLVER_H_ */

