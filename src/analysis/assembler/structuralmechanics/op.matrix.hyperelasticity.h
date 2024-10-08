
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
        SIMD C1 = load1(1.);
        SIMD F[9]; F[0] = C1; F[4] = C1; F[8] = C1;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (size_t n = 0; n < nodes; ++n) {
                    F[i * 3 + j] = F[i * 3 + j] + element.displacement[n][i] * element.dND[n][j];
                }
            }
        }
        printf("F:\n"); print(3, 3, F);
    }
};

}




#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_MATRIX_HYPERELASTICITY_H_ */
