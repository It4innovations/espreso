
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/math.h"
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
        SIMD JC[9];
        for (size_t n = 0; n < nodes; ++n) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    JC[i * 3 + j] = JC[i * 3 + j] + load1(element.dN[gp][n][i]) * (element.coords.node[n][j] + element.displacement[n][j]);
                }
            }
        }

        SIMD F[9]; multAtB<3>(F, JC, element.invJ, 1.0);
        SIMD eHat[9]; multAtB<3>(eHat, F, F, 0.5);
        SIMD eVec[6], C05 = load1(0.5), C2 = load1(2);
        eVec[0] = -C05 + eHat[0];
        eVec[1] = -C05 + eHat[4];
        eVec[2] = -C05 + eHat[8];
        eVec[3] =   C2 * eHat[1];
        eVec[4] =   C2 * eHat[5];
        eVec[5] =   C2 * eHat[2];
        SIMD sVec[6]; multMv<6, 6>(sVec, element.elasticity, eVec, load1(1.0));

        SIMD BL[6 * 3 * nodes];
        for (size_t n = 0; n < nodes; n++) {
            for (int j = 0; j < 3; j++) {
                BL[0 * 3 * nodes + n + j * nodes] = F[j * 3 + 0] * element.dND[n][0];
                BL[1 * 3 * nodes + n + j * nodes] = F[j * 3 + 1] * element.dND[n][1];
                BL[2 * 3 * nodes + n + j * nodes] = F[j * 3 + 2] * element.dND[n][2];
                BL[3 * 3 * nodes + n + j * nodes] = F[j * 3 + 0] * element.dND[n][1] + F[j * 3 + 1] * element.dND[n][0];
                BL[4 * 3 * nodes + n + j * nodes] = F[j * 3 + 1] * element.dND[n][2] + F[j * 3 + 2] * element.dND[n][1];
                BL[5 * 3 * nodes + n + j * nodes] = F[j * 3 + 0] * element.dND[n][2] + F[j * 3 + 2] * element.dND[n][0];
            }
        }

        SIMD scale = element.det * load1(element.w[gp]);
        // BLt * C * BL
        multAtBA<6, 3 * nodes>(element.K, BL, element.elasticity, scale);

        // GLt * S * GL -- block diagonal (TODO: create function)
        for (size_t n = 0; n < nodes; ++n) {
            SIMD a = element.dND[n][0] * sVec[0] + element.dND[n][1] * sVec[3] + element.dND[n][2] * sVec[5];
            SIMD b = element.dND[n][0] * sVec[3] + element.dND[n][1] * sVec[1] + element.dND[n][2] * sVec[4];
            SIMD c = element.dND[n][0] * sVec[5] + element.dND[n][1] * sVec[4] + element.dND[n][2] * sVec[2];

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

        // BLt * sVec
        multMtv<6, 3 * nodes>(element.nf, BL, sVec, scale);
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_ */
