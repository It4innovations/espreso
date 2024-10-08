
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_HYPERELASTICITY_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_HYPERELASTICITY_H_

#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/physics/structuralmechanics.h"

namespace espreso {

struct HyperElasticity: SubKernel {
    const char* name() const { return "Elasticity"; }

    HyperElasticity()
    : behaviour(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::PLANE_STRAIN),
      hyperElasticity(nullptr)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
    }

    void activate(StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour, const HyperElasticPropertiesConfiguration *hyperElasticity)
    {
        this->behaviour = behaviour;
        this->hyperElasticity = hyperElasticity;
        this->isconst = false;
        this->isactive = true;
    }

    StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR behaviour;
    const HyperElasticPropertiesConfiguration *hyperElasticity;
};

template <size_t ndim> struct HyperElasticityKernel;

template <> struct HyperElasticityKernel<2>: HyperElasticity {
    HyperElasticityKernel(const HyperElasticity &base): HyperElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (hyperElasticity->model) {
        case HyperElasticPropertiesConfiguration::MODEL::ARRUDA_BOYCE:
        case HyperElasticPropertiesConfiguration::MODEL::BLATZ_KO_FOAM:
        case HyperElasticPropertiesConfiguration::MODEL::GENT:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_2PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_3PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_5PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_9PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_CMP:
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_INC:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_1:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_2:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_3:
        default:
            break;
        }
    }
};

template <> struct HyperElasticityKernel<3>: HyperElasticity {
    HyperElasticityKernel(const HyperElasticity &base): HyperElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD vC[6], lambda;

        switch (hyperElasticity->model) {
        case HyperElasticPropertiesConfiguration::MODEL::ARRUDA_BOYCE:
        case HyperElasticPropertiesConfiguration::MODEL::BLATZ_KO_FOAM:
        case HyperElasticPropertiesConfiguration::MODEL::GENT:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_2PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_3PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_5PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_9PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_CMP:
        {
            // Voigt   0  1  2  3  4  5
            //        11 22 33 12 13 23
            SIMD C05 = load1(.5);
            SIMD J2 = vC[0] * vC[1] * vC[2] + load1(2) * vC[3] * vC[4] * vC[5] - vC[0] * vC[5] * vC[5] - vC[1] * vC[4] * vC[4] - vC[2] * vC[3] * vC[3];
//            SIMD J  = sqrt(J2);
            SIMD rJ2 = load1(1) / J2;

            SIMD vCinv[6] = {
                (vC[1] * vC[2] - vC[5] * vC[5] ) * rJ2,
                (vC[0] * vC[2] - vC[4] * vC[4] ) * rJ2,
                (vC[0] * vC[1] - vC[3] * vC[3] ) * rJ2,
                (vC[4] * vC[5] - vC[2] * vC[3] ) * rJ2,
                (vC[3] * vC[5] - vC[1] * vC[4] ) * rJ2,
                (vC[3] * vC[4] - vC[0] * vC[5] ) * rJ2,
            };
            SIMD symodot_vCinv[36] = {
                (vCinv[1] * vCinv[1] + vCinv[1] * vCinv[1]) * C05, (vCinv[4] * vCinv[4] + vCinv[4] * vCinv[4]) * C05, (vCinv[5] * vCinv[5] + vCinv[5] * vCinv[5]) * C05, (vCinv[1] * vCinv[4] + vCinv[4] * vCinv[1]) * C05, (vCinv[1] * vCinv[5] + vCinv[5] * vCinv[1]) * C05, (vCinv[4] * vCinv[5] + vCinv[5] * vCinv[4]) * C05,
                (vCinv[4] * vCinv[4] + vCinv[4] * vCinv[4]) * C05, (vCinv[2] * vCinv[2] + vCinv[2] * vCinv[2]) * C05, (vCinv[6] * vCinv[6] + vCinv[6] * vCinv[6]) * C05, (vCinv[4] * vCinv[2] + vCinv[2] * vCinv[4]) * C05, (vCinv[4] * vCinv[6] + vCinv[6] * vCinv[4]) * C05, (vCinv[2] * vCinv[6] + vCinv[6] * vCinv[2]) * C05,
                (vCinv[5] * vCinv[5] + vCinv[5] * vCinv[5]) * C05, (vCinv[6] * vCinv[6] + vCinv[6] * vCinv[6]) * C05, (vCinv[3] * vCinv[3] + vCinv[3] * vCinv[3]) * C05, (vCinv[5] * vCinv[6] + vCinv[6] * vCinv[5]) * C05, (vCinv[5] * vCinv[3] + vCinv[3] * vCinv[5]) * C05, (vCinv[6] * vCinv[3] + vCinv[3] * vCinv[6]) * C05,
                (vCinv[1] * vCinv[4] + vCinv[1] * vCinv[4]) * C05, (vCinv[4] * vCinv[2] + vCinv[4] * vCinv[2]) * C05, (vCinv[5] * vCinv[6] + vCinv[5] * vCinv[6]) * C05, (vCinv[1] * vCinv[2] + vCinv[4] * vCinv[4]) * C05, (vCinv[1] * vCinv[6] + vCinv[5] * vCinv[4]) * C05, (vCinv[4] * vCinv[6] + vCinv[5] * vCinv[2]) * C05,
                (vCinv[1] * vCinv[5] + vCinv[1] * vCinv[5]) * C05, (vCinv[4] * vCinv[6] + vCinv[4] * vCinv[6]) * C05, (vCinv[5] * vCinv[3] + vCinv[5] * vCinv[3]) * C05, (vCinv[1] * vCinv[6] + vCinv[4] * vCinv[5]) * C05, (vCinv[1] * vCinv[3] + vCinv[5] * vCinv[5]) * C05, (vCinv[4] * vCinv[3] + vCinv[5] * vCinv[6]) * C05,
                (vCinv[4] * vCinv[5] + vCinv[4] * vCinv[5]) * C05, (vCinv[2] * vCinv[6] + vCinv[2] * vCinv[6]) * C05, (vCinv[6] * vCinv[3] + vCinv[6] * vCinv[3]) * C05, (vCinv[4] * vCinv[6] + vCinv[2] * vCinv[5]) * C05, (vCinv[4] * vCinv[3] + vCinv[6] * vCinv[5]) * C05, (vCinv[2] * vCinv[3] + vCinv[6] * vCinv[6]) * C05
            };
//            double I1 = vC[0] + vC[1] + vC[3];
//            double I2 = vC[0] * vC[1] + vC[0] * vC[2] + vC[1] * vC[2] - vC[3] * vC[3] - vC[4] * vC[4] - vC[5] * vC[5];
            SIMD logJ2 = log(J2);
            SIMD al3 = C05 * (lambda * logJ2) - element.ecf.poissonRatio[0];
            element.S[0] = al3 * vCinv[0] + element.ecf.poissonRatio[0];
            element.S[1] = al3 * vCinv[1] + element.ecf.poissonRatio[0];
            element.S[2] = al3 * vCinv[2] + element.ecf.poissonRatio[0];
            element.S[3] = al3 * vCinv[3];
            element.S[4] = al3 * vCinv[4];
            element.S[5] = al3 * vCinv[5];

            SIMD be6 = lambda;
            SIMD be7 = load1(2) * element.ecf.poissonRatio[0] - lambda * logJ2;
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    element.elasticity[i * 6 + j] = be6 * vCinv[i] * vCinv[j] + be7 * symodot_vCinv[i * 6 + j];
                }
            }
        }
        break;
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_INC:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_1:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_2:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_3:
        default:
            break;
        }
    }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_HYPERELASTICITY_H_ */
