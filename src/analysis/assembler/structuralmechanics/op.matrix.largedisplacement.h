
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_

#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct MatrixLargeDisplacement: SubKernel {
    const char* name() const { return "MatrixLargeDisplacement"; }

    MatrixLargeDisplacement()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate()
    {
        this->isactive = 1;
    }
};

template <size_t nodes, size_t ndim> struct MatrixLargeDisplacementKernel;

template <size_t nodes>
struct MatrixLargeDisplacementKernel<nodes, 2>: MatrixLargeDisplacement {
    MatrixLargeDisplacementKernel(const MatrixLargeDisplacement &base): MatrixLargeDisplacement(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {

    }
};

template <size_t nodes>
struct MatrixLargeDisplacementKernel<nodes, 3>: MatrixLargeDisplacement {
    MatrixLargeDisplacementKernel(const MatrixLargeDisplacement &base): MatrixLargeDisplacement(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD scale = element.det * load1(element.w[gp]);

        SIMD BL[6 * 3 * nodes];
        for (size_t n = 0; n < nodes; n++) {
            for (int j = 0; j < 3; j++) {
                BL[0 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n][0];
                BL[1 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 1] * element.dND[n][1];
                BL[2 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 2] * element.dND[n][2];
                BL[3 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n][1] + element.F[j * 3 + 1] * element.dND[n][0];
                BL[4 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 1] * element.dND[n][2] + element.F[j * 3 + 2] * element.dND[n][1];
                BL[5 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n][2] + element.F[j * 3 + 2] * element.dND[n][0];
            }
        }

        // BLt * C * BL
        for (size_t n = 0; n < 3 * nodes; ++n) {
            SIMD a = zeros(), b = zeros(), c = zeros(), d = zeros(), e = zeros(), f = zeros();
            for (int i = 0; i < 6; ++i) {
                a = a + BL[i * 3 * nodes + n] * element.elasticity[i * 6 + 0];
                b = b + BL[i * 3 * nodes + n] * element.elasticity[i * 6 + 1];
                c = c + BL[i * 3 * nodes + n] * element.elasticity[i * 6 + 2];
                d = d + BL[i * 3 * nodes + n] * element.elasticity[i * 6 + 3];
                e = e + BL[i * 3 * nodes + n] * element.elasticity[i * 6 + 4];
                f = f + BL[i * 3 * nodes + n] * element.elasticity[i * 6 + 5];
            }
            SIMD nm = a * BL[0 * 3 * nodes + n] + b * BL[1 * 3 * nodes + n] + c * BL[2 * 3 * nodes + n] + d * BL[3 * 3 * nodes + n] + e * BL[4 * 3 * nodes + n] + f * BL[5 * 3 * nodes + n];
            element.K[n * 3 * nodes + n] = element.K[n * 3 * nodes + n] + scale * nm;
            for (size_t m = n + 1; m < 3 * nodes; ++m) {
                SIMD nm = a * BL[0 * 3 * nodes + m] + b * BL[1 * 3 * nodes + m] + c * BL[2 * 3 * nodes + m] + d * BL[3 * 3 * nodes + m] + e * BL[4 * 3 * nodes + m] + f * BL[5 * 3 * nodes + m];
                element.K[n * 3 * nodes + m] = element.K[n * 3 * nodes + m] + scale * nm;
                element.K[m * 3 * nodes + n] = element.K[m * 3 * nodes + n] + scale * nm;
            }
        }

        for (size_t n = 0; n < 3 * nodes; ++n) {
            SIMD nf =
                    element.sVec[0] * BL[0 * 3 * nodes + n] +
                    element.sVec[1] * BL[1 * 3 * nodes + n] +
                    element.sVec[2] * BL[2 * 3 * nodes + n] +
                    element.sVec[3] * BL[3 * 3 * nodes + n] +
                    element.sVec[4] * BL[4 * 3 * nodes + n] +
                    element.sVec[5] * BL[5 * 3 * nodes + n];
            element.nf[n] = element.nf[n] + scale * nf;
        }

        // GLt * S * GL -- block diagonal
        for (size_t n = 0; n < nodes; ++n) {
            SIMD a = element.dND[n][0] * element.sVec[0] + element.dND[n][1] * element.sVec[3] + element.dND[n][2] * element.sVec[5];
            SIMD b = element.dND[n][0] * element.sVec[3] + element.dND[n][1] * element.sVec[1] + element.dND[n][2] * element.sVec[4];
            SIMD c = element.dND[n][0] * element.sVec[5] + element.dND[n][1] * element.sVec[4] + element.dND[n][2] * element.sVec[2];

            SIMD nm = a * element.dND[n][0] + b * element.dND[n][1] + c * element.dND[n][2];
            element.K[(n + 0 * nodes) * 3 * nodes + (n + 0 * nodes)] = element.K[(n + 0 * nodes) * 3 * nodes + (n + 0 * nodes)] + scale * nm;
            element.K[(n + 1 * nodes) * 3 * nodes + (n + 1 * nodes)] = element.K[(n + 1 * nodes) * 3 * nodes + (n + 1 * nodes)] + scale * nm;
            element.K[(n + 2 * nodes) * 3 * nodes + (n + 2 * nodes)] = element.K[(n + 2 * nodes) * 3 * nodes + (n + 2 * nodes)] + scale * nm;
            for (size_t m = n + 1; m < nodes; ++m) {
                SIMD nm = a * element.dND[m][0] + b * element.dND[m][1] + c * element.dND[m][2];
                element.K[(n + 0 * nodes) * 3 * nodes + (m + 0 * nodes)] = element.K[(n + 0 * nodes) * 3 * nodes + (m + 0 * nodes)] + scale * nm;
                element.K[(n + 1 * nodes) * 3 * nodes + (m + 1 * nodes)] = element.K[(n + 1 * nodes) * 3 * nodes + (m + 1 * nodes)] + scale * nm;
                element.K[(n + 2 * nodes) * 3 * nodes + (m + 2 * nodes)] = element.K[(n + 2 * nodes) * 3 * nodes + (m + 2 * nodes)] + scale * nm;
                element.K[(m + 0 * nodes) * 3 * nodes + (n + 0 * nodes)] = element.K[(m + 0 * nodes) * 3 * nodes + (n + 0 * nodes)] + scale * nm;
                element.K[(m + 1 * nodes) * 3 * nodes + (n + 1 * nodes)] = element.K[(m + 1 * nodes) * 3 * nodes + (n + 1 * nodes)] + scale * nm;
                element.K[(m + 2 * nodes) * 3 * nodes + (n + 2 * nodes)] = element.K[(m + 2 * nodes) * 3 * nodes + (n + 2 * nodes)] + scale * nm;
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_ */
