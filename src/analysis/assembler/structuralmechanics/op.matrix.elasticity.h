
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_ELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_ELASTICITY_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/math.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct MatrixElasticity: SubKernel {
    const char* name() const { return "MatrixElasticity"; }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;

    MatrixElasticity()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN), nonlinear(true)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour, bool nonlinear)
    {
        this->behaviour = behaviour;
        this->isactive = 1;
        this->nonlinear = nonlinear;
    }

    bool nonlinear;
};

template <size_t nodes, size_t ndim> struct MatrixElasticityKernel;

template <size_t nodes>
struct MatrixElasticityKernel<nodes, 2>: MatrixElasticity {
    MatrixElasticityKernel(const MatrixElasticity &base): MatrixElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (behaviour) {
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS:
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
            // B = dX  0
            //      0 dY
            //     dY dX
            if (nonlinear) {
                SIMD BL[3 * 2 * nodes];
                for (size_t n = 0; n < nodes; n++) {
                    for (int j = 0; j < 2; j++) {
                        BL[0 * 2 * nodes + n + j * nodes] = element.F[j * 2 + 0] * element.dND[n * 2 + 0];
                        BL[1 * 2 * nodes + n + j * nodes] = element.F[j * 2 + 1] * element.dND[n * 2 + 1];
                        BL[2 * 2 * nodes + n + j * nodes] = element.F[j * 2 + 0] * element.dND[n * 2 + 1] + element.F[j * 2 + 1] * element.dND[n * 2 + 0];
                    }
                }
                SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]);
                multAtBA<3, 2 * nodes>(element.K, BL, element.vC4, scale);
                multAtB<2 * nodes, 3, 1>(element.nf, BL, element.vS, scale);

                SIMD S[4]; voigt3ToMatrix22(element.vS, S);
                SIMD KS[nodes * nodes];
                multABAt<nodes, 2>(KS, element.dND, S, scale);
                for (size_t n = 0; n < nodes; ++n) {
                    for (size_t m = 0; m < nodes; ++m) {
                        element.K[(n + 0 * nodes) * 2 * nodes + m + 0 * nodes] = element.K[(n + 0 * nodes) * 2 * nodes + m + 0 * nodes] + KS[n * nodes + m];
                        element.K[(n + 1 * nodes) * 2 * nodes + m + 1 * nodes] = element.K[(n + 1 * nodes) * 2 * nodes + m + 1 * nodes] + KS[n * nodes + m];
                    }
                }
            } else {
                SIMD B[3 * 2 * nodes];
                for (size_t n = 0; n < nodes; n++) {
                    B[0 * 2 * nodes + 0 * nodes + n] = element.dND[n * 2 + 0];
                    B[1 * 2 * nodes + 1 * nodes + n] = element.dND[n * 2 + 1];
                    B[2 * 2 * nodes + 0 * nodes + n] = element.dND[n * 2 + 1]; B[2 * 2 * nodes + 1 * nodes + n] = element.dND[n * 2 + 0];
                }
                multAtBA<3, 2 * nodes>(element.K, B, element.vC4, element.det * load1(element.w[gp]));
            }
            break;
        case StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
            // B = dX  0
            //      0 dY
            //      C  0
            //     dY dX
            if (nonlinear) {
                eslog::error("implement NONLINEAR AXISYMMETRIC assembler\n");
            } else {
                SIMD B[4 * 2 * nodes];
                for (size_t n = 0; n < nodes; n++) {
                    B[0 * 2 * nodes + 0 * nodes + n] = element.dND[n * 2 + 0];
                    B[1 * 2 * nodes + 1 * nodes + n] = element.dND[n * 2 + 1];
                    B[2 * 2 * nodes + 0 * nodes + n] = load1(element.N[gp][n]) / element.coords.gp[0];
                    B[3 * 2 * nodes + 0 * nodes + n] = element.dND[n * 2 + 1]; B[3 * 2 * nodes + 1 * nodes + n] = element.dND[n * 2 + 0];
                }
                multAtBA<4, 2 * nodes>(element.K, B, element.vC4, element.det * load1(element.w[gp]) * load1(2 * M_PI) * element.coords.gp[0]);
            }
            break;
        }
    }
};

template <size_t nodes>
struct MatrixElasticityKernel<nodes, 3>: MatrixElasticity {
    MatrixElasticityKernel(const MatrixElasticity &base): MatrixElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        // B = dX  0  0
        //      0 dY  0
        //      0  0 dZ
        //     dY dX  0
        //      0 dZ dY
        //     dZ  0 dX

        if (nonlinear) {
            SIMD BL[6 * 3 * nodes];
            for (size_t n = 0; n < nodes; n++) {
                for (int j = 0; j < 3; j++) {
                    BL[0 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n * 3 + 0];
                    BL[1 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 1] * element.dND[n * 3 + 1];
                    BL[2 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 2] * element.dND[n * 3 + 2];
                    BL[3 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n * 3 + 1] + element.F[j * 3 + 1] * element.dND[n * 3 + 0];
                    BL[4 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 1] * element.dND[n * 3 + 2] + element.F[j * 3 + 2] * element.dND[n * 3 + 1];
                    BL[5 * 3 * nodes + n + j * nodes] = element.F[j * 3 + 0] * element.dND[n * 3 + 2] + element.F[j * 3 + 2] * element.dND[n * 3 + 0];
                }
            }
            SIMD scale = element.det * load1(element.w[gp]);
            multAtBA<6, 3 * nodes>(element.K, BL, element.vC4, scale);
            multAtB<3 * nodes, 6, 1>(element.nf, BL, element.vS, scale);

            SIMD S[9]; voigt6ToMatrix33(element.vS, S);
            SIMD KS[nodes * nodes];
            multABAt<nodes, 3>(KS, element.dND, S, scale);
            for (size_t n = 0; n < nodes; ++n) {
                for (size_t m = 0; m < nodes; ++m) {
                    element.K[(n + 0 * nodes) * 3 * nodes + m + 0 * nodes] = element.K[(n + 0 * nodes) * 3 * nodes + m + 0 * nodes] + KS[n * nodes + m];
                    element.K[(n + 1 * nodes) * 3 * nodes + m + 1 * nodes] = element.K[(n + 1 * nodes) * 3 * nodes + m + 1 * nodes] + KS[n * nodes + m];
                    element.K[(n + 2 * nodes) * 3 * nodes + m + 2 * nodes] = element.K[(n + 2 * nodes) * 3 * nodes + m + 2 * nodes] + KS[n * nodes + m];
                }
            }
        } else {
            SIMD B[6 * 3 * nodes];
            for (size_t n = 0; n < nodes; n++) {
                B[0 * 3 * nodes + 0 * nodes + n] = element.dND[n * 3 + 0];
                B[1 * 3 * nodes + 1 * nodes + n] = element.dND[n * 3 + 1];
                B[2 * 3 * nodes + 2 * nodes + n] = element.dND[n * 3 + 2];
                B[3 * 3 * nodes + 0 * nodes + n] = element.dND[n * 3 + 1]; B[3 * 3 * nodes + 1 * nodes + n] = element.dND[n * 3 + 0];
                B[4 * 3 * nodes + 1 * nodes + n] = element.dND[n * 3 + 2]; B[4 * 3 * nodes + 2 * nodes + n] = element.dND[n * 3 + 1];
                B[5 * 3 * nodes + 0 * nodes + n] = element.dND[n * 3 + 2]; B[5 * 3 * nodes + 2 * nodes + n] = element.dND[n * 3 + 0];
            }
            multAtBA<6, 3 * nodes>(element.K, B, element.vC4, element.det * load1(element.w[gp]));
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_ELASTICITY_H_ */
