
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_PRINT_EIGENVALUES_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_PRINT_EIGENVALUES_H_

#include "subkernel.h"
#include "math/primitives/matrix_dense.h"
#include "math/primitives/vector_dense.h"
#include "wrappers/simd/simd.h"

namespace espreso {

struct PrintEigenValues: SubKernel {
    const char* name() const { return "PrintEigenValues"; }

    PrintEigenValues()
    {
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate()
    {
        this->isactive = true;
    }
};

struct PrintEigenValuesKernel: PrintEigenValues {
    PrintEigenValuesKernel(const PrintEigenValues &base): PrintEigenValues(base) {}

    Matrix_Dense<double> m;
    Vector_Dense<double> e;

    void simd(const char * name, size_t size, SIMD M[], size_t n);
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_PRINT_EIGENVALUES_H_ */
