
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

template <size_t nodes, size_t ndim> struct HyperElasticityKernel;

template <size_t nodes> struct HyperElasticityKernel<nodes, 2>: HyperElasticity {
    HyperElasticityKernel(const HyperElasticity &base): HyperElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (hyperElasticity->model) {
        case HyperElasticPropertiesConfiguration::MODEL::KIRCHHOFF:
        case HyperElasticPropertiesConfiguration::MODEL::ARRUDA_BOYCE:
        case HyperElasticPropertiesConfiguration::MODEL::BLATZ_KO_FOAM:
        case HyperElasticPropertiesConfiguration::MODEL::GENT:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_2PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_3PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_5PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_9PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEAN_CMP:
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEAN_INC:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_1:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_2:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_3:
        default:
            break;
        }
    }
};

template <size_t nodes> struct HyperElasticityKernel<nodes, 3>: HyperElasticity {
    HyperElasticityKernel(const HyperElasticity &base): HyperElasticity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                element.F[i * 3 + j] = zeros();
                for (size_t n = 0; n < nodes; ++n) {
                    element.F[i * 3 + j] = element.F[i * 3 + j] + element.displacement[n][i] * element.dND[n * 3 + j];
                }
            }
        }
        element.F[0] = element.F[0] + load1(1.);
        element.F[4] = element.F[4] + load1(1.);
        element.F[8] = element.F[8] + load1(1.);

        // Voigt   0  1  2  3  4  5
        //        11 22 33 12 13 23
        SIMD C2[9]; multAtB<3, 3, 3>(C2, element.F, element.F); // Right Cauchy-Green tensor

        switch (hyperElasticity->model) {
        case HyperElasticPropertiesConfiguration::MODEL::KIRCHHOFF:
        {
            SIMD E = element.ecf.youngModulus[0];
            SIMD nu = element.ecf.poissonRatio[0]; // https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
            SIMD lambda = E * nu / ((load1(1.) + nu) * (load1(1.) - load1(2.) * nu));
            SIMD mu = E / (load1(2.) + load1(2.) * nu);
            element.elasticity[ 0] = element.elasticity[ 7] = element.elasticity[14] = lambda + load1(2.) * mu;
            element.elasticity[ 1] = element.elasticity[ 2] = element.elasticity[ 8] = lambda;
            element.elasticity[21] = element.elasticity[28] = element.elasticity[35] = mu;

            SIMD al0 = load1(.5) * (C2[0] + C2[4] + C2[8]) * lambda - mu - load1(.5) * load1(3) * lambda;
            element.vS[0] = al0 + mu * C2[0];
            element.vS[1] = al0 + mu * C2[4];
            element.vS[2] = al0 + mu * C2[8];
            element.vS[3] = mu * C2[1];
            element.vS[4] = mu * C2[5];
            element.vS[5] = mu * C2[2];
        }
        break;
        case HyperElasticPropertiesConfiguration::MODEL::ARRUDA_BOYCE:
        case HyperElasticPropertiesConfiguration::MODEL::BLATZ_KO_FOAM:
        case HyperElasticPropertiesConfiguration::MODEL::GENT:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_2PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_3PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_5PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_9PARAMS:
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEAN_CMP:
        {
            SIMD E = element.ecf.youngModulus[0];
            SIMD nu = element.ecf.poissonRatio[0]; // https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
            SIMD lambda = E * nu / ((load1(1.) + nu) * (load1(1.) - load1(2.) * nu));
            SIMD mu = E / (load1(2.) + load1(2.) * nu);

            SIMD J2 = C2[0] * C2[4] * C2[8] + load1(2) * C2[1] * C2[5] * C2[2] - C2[0] * C2[2] * C2[2] - C2[4] * C2[5] * C2[5] - C2[8] * C2[1] * C2[1];
            SIMD rJ2 = load1(1) / J2;
            SIMD vCinv[6] = {
                (C2[4] * C2[8] - C2[2] * C2[2] ) * rJ2,
                (C2[0] * C2[8] - C2[5] * C2[5] ) * rJ2,
                (C2[0] * C2[4] - C2[1] * C2[1] ) * rJ2,
                (C2[5] * C2[2] - C2[8] * C2[1] ) * rJ2,
                (C2[1] * C2[2] - C2[4] * C2[5] ) * rJ2,
                (C2[1] * C2[5] - C2[0] * C2[2] ) * rJ2,
            };

            SIMD C05 = load1(.5);
            SIMD logJ2 = log(J2);
            SIMD al3 = C05 * (lambda * logJ2) - mu;
            element.vS[0] = al3 * vCinv[0] + mu;
            element.vS[1] = al3 * vCinv[1] + mu;
            element.vS[2] = al3 * vCinv[2] + mu;
            element.vS[3] = al3 * vCinv[3];
            element.vS[4] = al3 * vCinv[4];
            element.vS[5] = al3 * vCinv[5];

            auto ij = [] (int i, int j) { return (i - 1) * 6 + j - 1; };
            SIMD symodot_vCinv[36];
            symodot_vCinv[ij(1,1)] =  vCinv[0] * vCinv[0];
            symodot_vCinv[ij(1,2)] =  vCinv[3] * vCinv[3];                                      symodot_vCinv[ij(2,1)] = symodot_vCinv[ij(1,2)];
            symodot_vCinv[ij(1,3)] =  vCinv[5] * vCinv[5];                                      symodot_vCinv[ij(3,1)] = symodot_vCinv[ij(1,3)];
            symodot_vCinv[ij(1,4)] =  vCinv[0] * vCinv[3];                                      symodot_vCinv[ij(4,1)] = symodot_vCinv[ij(1,4)];
            symodot_vCinv[ij(1,5)] =  vCinv[3] * vCinv[5];                                      symodot_vCinv[ij(5,1)] = symodot_vCinv[ij(1,5)];
            symodot_vCinv[ij(1,6)] =  vCinv[0] * vCinv[5];                                      symodot_vCinv[ij(6,1)] = symodot_vCinv[ij(1,6)];
            symodot_vCinv[ij(2,2)] =  vCinv[1] * vCinv[1];
            symodot_vCinv[ij(2,3)] =  vCinv[4] * vCinv[4];                                      symodot_vCinv[ij(3,2)] = symodot_vCinv[ij(2,3)];
            symodot_vCinv[ij(2,4)] =  vCinv[1] * vCinv[3];                                      symodot_vCinv[ij(4,2)] = symodot_vCinv[ij(2,4)];
            symodot_vCinv[ij(2,5)] =  vCinv[1] * vCinv[4];                                      symodot_vCinv[ij(5,2)] = symodot_vCinv[ij(2,5)];
            symodot_vCinv[ij(2,6)] =  vCinv[3] * vCinv[4];                                      symodot_vCinv[ij(6,2)] = symodot_vCinv[ij(2,6)];
            symodot_vCinv[ij(3,3)] =  vCinv[2] * vCinv[2];
            symodot_vCinv[ij(3,4)] =  vCinv[4] * vCinv[5];                                      symodot_vCinv[ij(4,3)] = symodot_vCinv[ij(3,4)];
            symodot_vCinv[ij(3,5)] =  vCinv[2] * vCinv[4];                                      symodot_vCinv[ij(5,3)] = symodot_vCinv[ij(3,5)];
            symodot_vCinv[ij(3,6)] =  vCinv[2] * vCinv[5];                                      symodot_vCinv[ij(6,3)] = symodot_vCinv[ij(3,6)];
            symodot_vCinv[ij(4,4)] = (vCinv[0] * vCinv[1] + vCinv[3] * vCinv[3]) * load1(.5);
            symodot_vCinv[ij(4,5)] = (vCinv[3] * vCinv[4] + vCinv[1] * vCinv[5]) * load1(.5);   symodot_vCinv[ij(5,4)] = symodot_vCinv[ij(4,5)];
            symodot_vCinv[ij(4,6)] = (vCinv[0] * vCinv[4] + vCinv[3] * vCinv[5]) * load1(.5);   symodot_vCinv[ij(6,4)] = symodot_vCinv[ij(4,6)];
            symodot_vCinv[ij(5,5)] = (vCinv[1] * vCinv[2] + vCinv[4] * vCinv[4]) * load1(.5);
            symodot_vCinv[ij(5,6)] = (vCinv[2] * vCinv[3] + vCinv[4] * vCinv[5]) * load1(.5);   symodot_vCinv[ij(6,5)] = symodot_vCinv[ij(5,6)];
            symodot_vCinv[ij(6,6)] = (vCinv[0] * vCinv[2] + vCinv[5] * vCinv[5]) * load1(.5);

            SIMD be6 = lambda;
            SIMD be7 = load1(2) * mu - lambda * logJ2;
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    element.elasticity[i * 6 + j] = be6 * vCinv[i] * vCinv[j] + be7 * symodot_vCinv[i * 6 + j];
                }
            }
        }
        break;
        case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEAN_INC:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_1:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_2:
        case HyperElasticPropertiesConfiguration::MODEL::OGDEN_3:
        default:
            break;
        }

        // make the elasticity symmetric
        element.elasticity[ 6] = element.elasticity[ 1];
        element.elasticity[12] = element.elasticity[ 2]; element.elasticity[13] = element.elasticity[ 8];
        element.elasticity[18] = element.elasticity[ 3]; element.elasticity[19] = element.elasticity[ 9]; element.elasticity[20] = element.elasticity[15];
        element.elasticity[24] = element.elasticity[ 4]; element.elasticity[25] = element.elasticity[10]; element.elasticity[26] = element.elasticity[16]; element.elasticity[27] = element.elasticity[22];
        element.elasticity[30] = element.elasticity[ 5]; element.elasticity[31] = element.elasticity[11]; element.elasticity[32] = element.elasticity[17]; element.elasticity[33] = element.elasticity[23]; element.elasticity[34] = element.elasticity[29];
    }
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_HYPERELASTICITY_H_ */
