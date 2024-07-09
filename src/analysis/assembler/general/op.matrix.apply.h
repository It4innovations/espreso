
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_MATRIX_APPLY_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_MATRIX_APPLY_H_

#include "element.h"
#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct MatrixApply: SubKernel {
    const char* name() const { return "MatrixApply"; }

    MatrixApply()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate()
    {
        this->isactive = 1;
    }
};

template <size_t dimension>
struct MatrixApplyKernel: MatrixApply {
    MatrixApplyKernel(const MatrixApply &base): MatrixApply(base) {}

    void simd(SIMD out[dimension], SIMD matrix[dimension * dimension], SIMD in[dimension], const double &scale)
    {
        SIMD s = load1(scale);
        for (size_t i = 0; i < dimension; ++i) {
            for (size_t j = 0; j < dimension; ++j) {
                out[i] = out[i] + s * matrix[i * dimension + j] * in[j];
            }
        }
    }
};

}


#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_MATRIX_APPLY_H_ */
