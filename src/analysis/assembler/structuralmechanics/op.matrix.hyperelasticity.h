
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_HYPERELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_HYPERELASTICITY_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/math.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct MatrixHyperElasticity: SubKernel {
    const char* name() const { return "MatrixHyperElasticity"; }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;

    MatrixHyperElasticity()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour)
    {
        this->behaviour = behaviour;
        this->isactive = 1;
    }
};

template <size_t nodes, size_t ndim> struct MatrixHyperElasticityKernel;

template <size_t nodes>
struct MatrixHyperElasticityKernel<nodes, 2>: MatrixHyperElasticity {
    MatrixHyperElasticityKernel(const MatrixHyperElasticity &base): MatrixHyperElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {

    }
};

template <size_t nodes>
struct MatrixHyperElasticityKernel<nodes, 3>: MatrixHyperElasticity {
    MatrixHyperElasticityKernel(const MatrixHyperElasticity &base): MatrixHyperElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
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
        multAtBA<6, 3 * nodes>(element.K, BL, element.elasticity, scale);
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
    }
};

}




#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_HYPERELASTICITY_H_ */
