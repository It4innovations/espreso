
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_MASS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_MASS_H_

#include "analysis/assembler/general/subkernel.h"
#include "math/primitives/matrix_info.h"
#include "wrappers/simd/simd.h"

namespace espreso {

struct MatrixMass: public SubKernel {
    const char* name() const { return "MatrixMass"; }

    MatrixMass()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate()
    {
        this->isactive = 1;
    }
};

template <size_t nodes, size_t dim> struct MatrixMassKernel;

template <size_t nodes>
struct MatrixMassKernel<nodes, 1>: MatrixMass {
    MatrixMassKernel(const MatrixMass &base): MatrixMass(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD scale = element.det * load1(element.w[gp]) * element.ecf.density;
        for (size_t n1 = 0; n1 < nodes; ++n1) {
            element.M[n1 * nodes + n1] = element.M[n1 * nodes + n1] + scale * load1(element.N[gp][n1] * element.N[gp][n1]);
            for (size_t m1 = n1 + 1; m1 < nodes; ++m1) {
                SIMD k = scale * load1(element.N[gp][n1] * element.N[gp][m1]);
                element.M[n1 * nodes + m1] = element.M[n1 * nodes + m1] + k;
                element.M[m1 * nodes + n1] = element.M[m1 * nodes + n1] + k;
            }
        }
    }
};

template <size_t nodes>
struct MatrixMassKernel<nodes, 2>: MatrixMass {
    MatrixMassKernel(const MatrixMass &base): MatrixMass(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD scale = element.det * load1(element.w[gp]) * element.ecf.density;
        for (size_t n1 = 0, n2 = nodes; n1 < nodes; ++n1, ++n2) {
            element.M[n1 * 2 * nodes + n1] = element.M[n1 * 2 * nodes + n1] + scale * load1(element.N[gp][n1] * element.N[gp][n1]);
            element.M[n2 * 2 * nodes + n2] = element.M[n2 * 2 * nodes + n2] + scale * load1(element.N[gp][n1] * element.N[gp][n1]);
            for (size_t m1 = n1 + 1, m2 = n2 + 1; m1 < nodes; ++m1, ++m2) {
                SIMD k = scale * load1(element.N[gp][n1] * element.N[gp][m1]);
                element.M[n1 * 2 * nodes + m1] = element.M[n1 * 2 * nodes + m1] + k;
                element.M[n2 * 2 * nodes + m2] = element.M[n2 * 2 * nodes + m2] + k;
                element.M[m1 * 2 * nodes + n1] = element.M[m1 * 2 * nodes + n1] + k;
                element.M[m2 * 2 * nodes + n2] = element.M[m2 * 2 * nodes + n2] + k;
            }
        }
    }
};

template <size_t nodes>
struct MatrixMassKernel<nodes, 3>: MatrixMass {
    MatrixMassKernel(const MatrixMass &base): MatrixMass(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD scale = element.det * load1(element.w[gp]) * element.ecf.density;
        for (size_t n1 = 0, n2 = nodes, n3 = 2 * nodes; n1 < nodes; ++n1, ++n2, ++n3) {
            element.M[n1 * 3 * nodes + n1] = element.M[n1 * 3 * nodes + n1] + scale * load1(element.N[gp][n1] * element.N[gp][n1]);
            element.M[n2 * 3 * nodes + n2] = element.M[n2 * 3 * nodes + n2] + scale * load1(element.N[gp][n1] * element.N[gp][n1]);
            element.M[n3 * 3 * nodes + n3] = element.M[n3 * 3 * nodes + n3] + scale * load1(element.N[gp][n1] * element.N[gp][n1]);
            for (size_t m1 = n1 + 1, m2 = n2 + 1, m3 = n3 + 1; m1 < nodes; ++m1, ++m2, ++m3) {
                SIMD k = scale * load1(element.N[gp][n1] * element.N[gp][m1]);
                element.M[n1 * 3 * nodes + m1] = element.M[n1 * 3 * nodes + m1] + k;
                element.M[n2 * 3 * nodes + m2] = element.M[n2 * 3 * nodes + m2] + k;
                element.M[n3 * 3 * nodes + m3] = element.M[n3 * 3 * nodes + m3] + k;
                element.M[m1 * 3 * nodes + n1] = element.M[m1 * 3 * nodes + n1] + k;
                element.M[m2 * 3 * nodes + n2] = element.M[m2 * 3 * nodes + n2] + k;
                element.M[m3 * 3 * nodes + n3] = element.M[m3 * 3 * nodes + n3] + k;
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_MASS_H_ */
