
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/math.h"
#include "config/ecf/physics/structuralmechanics.h"


namespace espreso {

struct MatrixLargeDisplacement: SubKernel {
    const char* name() const { return "MatrixLargeDisplacement"; }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
    LinearElasticPropertiesConfiguration::MODEL model;

    MatrixLargeDisplacement()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN),
      model(LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour, LinearElasticPropertiesConfiguration::MODEL model)
    {
        this->behaviour = behaviour;
        this->model = model;
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
        for (size_t n = 0; n < nodes; ++n) {

        }

        SIMD JC[4];
        for (size_t n = 0; n < nodes; ++n) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    JC[i * 2 + j] = JC[i * 2 + j] + load1(element.dN[gp][n][i]) * (element.coords.node[n][j] + element.displacement[n][j]);
                }
            }
        }

        SIMD F[4]; multAtB<2, 2, 2>(F, JC, element.invJ, load1(1.0));
        SIMD eHat[4]; multAtB<2, 2, 2>(eHat, F, F, load1(.5));
        SIMD eVec[3], C05 = load1(0.5);
        eVec[0] = -C05 + eHat[0];
        eVec[1] = -C05 + eHat[3];
        eVec[2] = eHat[1] + eHat[2];
        SIMD sVec[3];
        multAB<3, 3, 1>(sVec, element.elasticity, eVec, load1(1.0));

//        SIMD fbar[12] = {
//            F(1,1) , zeros(), F(2,1) , zeros(),
//            zeros(), F(1,2) , zeros(), F(2,2) ,
//            F(1,2) , F(1,1) , F(2,2) , F(2,1) ,
//        };
        SIMD fbar[12] = {
            F[0]   , zeros(), F[2]   , zeros(),
            zeros(), F[1]   , zeros(), F[3]   ,
            F[1]   , F[0]   , F[3]   , F[2]   ,
        };
//        SIMD bt[4 * 2 * nodes] = {
//            dnx(1), 0     , dnx(2), 0     , dnx(3), 0     , dnx(4), 0     ,
//            dny(1), 0     , dny(2), 0     , dny(3), 0     , dny(4), 0     ,
//            0     , dnx(1), 0     , dnx(2), 0     , dnx(3), 0     , dnx(4),
//            0     , dny(1), 0     , dny(2), 0     , dny(3), 0     , dny(4),
//        };
        SIMD bt[4 * 2 * nodes];
        for (size_t n = 0; n < nodes; ++n) {
            bt[0 * 2 * nodes + 2 * n] = element.dND[n][0];
            bt[1 * 2 * nodes + 2 * n] = element.dND[n][1];
            bt[2 * 2 * nodes + 2 * n + 1] = element.dND[n][0];
            bt[3 * 2 * nodes + 2 * n + 1] = element.dND[n][1];
        }

        // GT = fbar * bt
        SIMD GT[3 * 2 * nodes]; multAB<3, 4, 2 * nodes>(GT, fbar, bt);
//        for (size_t n = 0; n < nodes; n++) {
//            GT[0 * 2 * nodes + n        ] = F[0] * element.dND[n][0];
//            GT[0 * 2 * nodes + n + nodes] = F[2] * element.dND[n][0];
//            GT[1 * 2 * nodes + n        ] = F[1] * element.dND[n][1];
//            GT[1 * 2 * nodes + n + nodes] = F[3] * element.dND[n][1];
//            GT[2 * 2 * nodes + n        ] = F[1] * element.dND[n][0] + F[0] * element.dND[n][1];
//            GT[2 * 2 * nodes + n + nodes] = F[3] * element.dND[n][0] + F[1] * element.dND[n][1];
//        }

        SIMD KK[2 * nodes * 2 * nodes];

        SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]);
        // G * C * GT
        multAtBA<3, 2 * nodes>(KK, GT, element.elasticity, scale);

        // b * sbar * bt -- block diagonal (TODO: create function)
        SIMD sbar[16];
        sbar[0] = sVec[0];
        sbar[1] = sVec[2];
        sbar[4] = sVec[2];
        sbar[5] = sVec[1];

        sbar[10] = sVec[0];
        sbar[11] = sVec[2];
        sbar[14] = sVec[2];
        sbar[15] = sVec[1];

        multAtBA<4, 2 * nodes>(KK, bt, sbar, scale);

        for (size_t r = 0; r < nodes; ++r) {
            for (size_t c = 0; c < nodes; ++c) {
                element.K[r * 2 * nodes                     + c]         = element.K[r * 2 * nodes                     + c]         + KK[(2 * r)     * 2 * nodes + (2 * c)];
                element.K[r * 2 * nodes                     + c + nodes] = element.K[r * 2 * nodes                     + c + nodes] + KK[(2 * r)     * 2 * nodes + (2 * c + 1)];
                element.K[r * 2 * nodes + nodes * 2 * nodes + c]         = element.K[r * 2 * nodes + nodes * 2 * nodes + c]         + KK[(2 * r + 1) * 2 * nodes + (2 * c)];
                element.K[r * 2 * nodes + nodes * 2 * nodes + c + nodes] = element.K[r * 2 * nodes + nodes * 2 * nodes + c + nodes] + KK[(2 * r + 1) * 2 * nodes + (2 * c + 1)];
            }
        }

//        for (size_t n = 0; n < nodes; ++n) {
//            SIMD a = element.dND[n][0] * sVec[0] + element.dND[n][1] * sVec[2];
//            SIMD b = element.dND[n][0] * sVec[2] + element.dND[n][1] * sVec[1];
//
//            SIMD nm = a * element.dND[n][0] + b * element.dND[n][1];
//            element.K[(n + 0 * nodes) * 2 * nodes + (n + 0 * nodes)] = element.K[(n + 0 * nodes) * 2 * nodes + (n + 0 * nodes)] + scale * nm;
//            element.K[(n + 1 * nodes) * 2 * nodes + (n + 1 * nodes)] = element.K[(n + 1 * nodes) * 2 * nodes + (n + 1 * nodes)] + scale * nm;
//            for (size_t m = n + 1; m < nodes; ++m) {
//                SIMD nm = a * element.dND[m][0] + b * element.dND[m][1];
//                element.K[(n + 0 * nodes) * 2 * nodes + (m + 0 * nodes)] = element.K[(n + 0 * nodes) * 2 * nodes + (m + 0 * nodes)] + scale * nm;
//                element.K[(n + 1 * nodes) * 2 * nodes + (m + 1 * nodes)] = element.K[(n + 1 * nodes) * 2 * nodes + (m + 1 * nodes)] + scale * nm;
//                element.K[(m + 0 * nodes) * 2 * nodes + (n + 0 * nodes)] = element.K[(m + 0 * nodes) * 2 * nodes + (n + 0 * nodes)] + scale * nm;
//                element.K[(m + 1 * nodes) * 2 * nodes + (n + 1 * nodes)] = element.K[(m + 1 * nodes) * 2 * nodes + (n + 1 * nodes)] + scale * nm;
//            }
//        }

        SIMD RR[2 * nodes];
        // G * sVec
        multAtB<3, 2 * nodes, 1>(RR, GT, sVec, scale);

        for (size_t r = 0; r < nodes; ++r) {
            element.nf[r]         = element.nf[r]         + RR[2 * r];
            element.nf[r + nodes] = element.nf[r + nodes] + RR[2 * r + 1];
        }
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

        SIMD F[9]; multAtB<3, 3, 3>(F, JC, element.invJ, load1(1.0));
        SIMD eHat[9]; multAtB<3, 3, 3>(eHat, F, F, load1(.5));
        SIMD eVec[6], C05 = load1(0.5), C2 = load1(2);
        eVec[0] = -C05 + eHat[0];
        eVec[1] = -C05 + eHat[4];
        eVec[2] = -C05 + eHat[8];
        eVec[3] =   C2 * eHat[1];
        eVec[4] =   C2 * eHat[5];
        eVec[5] =   C2 * eHat[2];
        SIMD sVec[6];
        multAB<6, 6, 1>(sVec, element.elasticity, eVec, load1(1.0));

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
        multAtB<6, 3 * nodes, 1>(element.nf, BL, sVec, scale);
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_MATRIX_LARGEDISPLACEMENT_H_ */
