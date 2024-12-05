
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_

#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/math.h"
#include "analysis/assembler/general/op.integration.h"
#include "analysis/assembler/structuralmechanics/op.material.h"
#include "basis/containers/serializededata.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nameddata.h"

namespace espreso {

struct Stress: SubKernel {
    const char* name() const { return "Stress"; }

    Stress();
    serializededata<esint, esint>::const_iterator enodes, end;
    double *multiplicity;
    double *principalStress, *principalStrain, *componentStress, *componentStrain, *vonMisesStress, *vonMisesStrain;
    double *principalStressAvg, *principalStrainAvg, *componentStressAvg, *componentStrainAvg, *vonMisesStressAvg, *vonMisesStrainAvg;

    void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, size_t interval, double *multiplicity);
};

template <size_t nodes, size_t gps, size_t ndim> struct StressKernel;

template <size_t nodes, size_t gps>
struct StressKernel<nodes, gps, 2>: Stress {
    StressKernel(const Stress &base): Stress(base) {}

    template <typename Element>
    void init(Element &element)
    {
        for (int e = 0; e < element.elements; ++e) {
            principalStress[2 * e + 0] = 0;
            principalStress[2 * e + 1] = 0;
            principalStrain[2 * e + 0] = 0;
            principalStrain[2 * e + 1] = 0;

            componentStress[3 * e + 0] = 0;
            componentStress[3 * e + 1] = 0;
            componentStress[3 * e + 2] = 0;
            componentStrain[3 * e + 0] = 0;
            componentStrain[3 * e + 1] = 0;
            componentStrain[3 * e + 2] = 0;

            vonMisesStress[e] = 0;
            vonMisesStrain[e] = 0;
        }
    }

    template <typename Element>
    void move(Element &element)
    {
        principalStress += SIMD::size * 2;
        principalStrain += SIMD::size * 2;
        componentStress += SIMD::size * 3;
        componentStrain += SIMD::size * 3;
        vonMisesStress  += SIMD::size;
        vonMisesStrain  += SIMD::size;

        for (int e = 0; e < element.elements; ++e) {
            ++enodes;
        }
    }

    template <typename Element>
    void element(Element &element)
    {
        constexpr double rgps = 1.0 / gps;
        SIMD eig[2], nu = element.ecf.poissonRatio[0];
        SIMD S[4] = {
                element.vS[0], element.vS[2],
                element.vS[2], element.vS[1],
        };
        eigSym22Desc(S, eig);
        for (int e = 0; e < element.elements; ++e) {
            principalStress[2 * e + 0] += rgps * eig[0][e];
            principalStress[2 * e + 1] += rgps * eig[1][e];

            componentStress[3 * e + 0] += rgps * S[0][e];
            componentStress[3 * e + 1] += rgps * S[3][e];
            componentStress[3 * e + 2] += rgps * S[1][e];

            vonMisesStress[e] += rgps * sqrt(.5 * ((eig[0][e] - eig[1][e]) * (eig[0][e] - eig[1][e])));
        }

        SIMD E[4] = {
                load1(.5) * element.C2[0] - load1(.5), load1(.5) * element.C2[1]            ,
                load1(.5) * element.C2[2]            , load1(.5) * element.C2[3] - load1(.5),
        };
        eigSym22Desc(E, eig);
        for (int e = 0; e < element.elements; ++e) {
            principalStrain[2 * e + 0] += rgps * eig[0][e];
            principalStrain[2 * e + 1] += rgps * eig[1][e];

            componentStrain[3 * e + 0] += rgps * E[0][e];
            componentStrain[3 * e + 1] += rgps * E[3][e];
            componentStrain[3 * e + 2] += rgps * E[1][e];

            vonMisesStrain[e] += rgps * sqrt(.5 * ((eig[0][e] - eig[1][e]) * (eig[0][e] - eig[1][e]))) / (1 + nu[e]);
        }
    }

    template <typename Element>
    void node(Element &element, size_t n)
    {
        size_t index = 0;
        for (size_t nn = 0; nn < nodes; ++nn) {
            index += nn * element.NN[n][nn]; // get node index, only one values in NN is 1, other are 0
        }
        SIMD nu = element.ecf.poissonRatio[0];
        SIMD eigS[2], S[4] = {
                element.vS[0], element.vS[2],
                element.vS[2], element.vS[1],
        };
        eigSym22Desc(S, eigS);
        SIMD eigE[2], E[4] = {
                load1(.5) * element.C2[0] - load1(.5), load1(.5) * element.C2[1]            ,
                load1(.5) * element.C2[2]            , load1(.5) * element.C2[3] - load1(.5),
        };
        eigSym22Desc(E, eigE);
        auto enodes = this->enodes;
        for (int e = 0; e < element.elements; ++e, ++enodes) {
            size_t node = enodes->at(index);
            double scale = multiplicity[node];
            principalStressAvg[2 * node + 0] += scale * eigS[0][e];
            principalStressAvg[2 * node + 1] += scale * eigS[1][e];
            componentStressAvg[3 * node + 0] += scale * S[0][e];
            componentStressAvg[3 * node + 1] += scale * S[3][e];
            componentStressAvg[3 * node + 2] += scale * S[1][e];
            vonMisesStressAvg [    node    ] += scale * sqrt(.5 * ((eigS[0][e] - eigS[1][e]) * (eigS[0][e] - eigS[1][e])));

            principalStrainAvg[2 * node + 0] += scale * eigE[0][e];
            principalStrainAvg[2 * node + 1] += scale * eigE[1][e];
            componentStrainAvg[3 * node + 0] += scale * E[0][e];
            componentStrainAvg[3 * node + 1] += scale * E[3][e];
            componentStrainAvg[3 * node + 2] += scale * E[1][e];
            vonMisesStrainAvg [    node    ] += scale * sqrt(.5 * ((eigE[0][e] - eigE[1][e]) * (eigE[0][e] - eigE[1][e]))) / (1 + nu[e]);
        }
    }
};

template <size_t nodes, size_t gps>
struct StressKernel<nodes, gps, 3>: Stress {
    StressKernel(const Stress &base): Stress(base) {}

    template <typename Element>
    void init(Element &element)
    {
        for (int e = 0; e < element.elements; ++e) {
            principalStress[3 * e + 0] = 0;
            principalStress[3 * e + 1] = 0;
            principalStress[3 * e + 2] = 0;
            principalStrain[3 * e + 0] = 0;
            principalStrain[3 * e + 1] = 0;
            principalStrain[3 * e + 2] = 0;

            componentStress[6 * e + 0] = 0;
            componentStress[6 * e + 1] = 0;
            componentStress[6 * e + 2] = 0;
            componentStress[6 * e + 3] = 0;
            componentStress[6 * e + 4] = 0;
            componentStress[6 * e + 5] = 0;
            componentStrain[6 * e + 0] = 0;
            componentStrain[6 * e + 1] = 0;
            componentStrain[6 * e + 2] = 0;
            componentStrain[6 * e + 3] = 0;
            componentStrain[6 * e + 4] = 0;
            componentStrain[6 * e + 5] = 0;

            vonMisesStress[e] = 0;
            vonMisesStrain[e] = 0;
        }
    }

    template <typename Element>
    void move(Element &element)
    {
        principalStress += SIMD::size * 3;
        principalStrain += SIMD::size * 3;
        componentStress += SIMD::size * 6;
        componentStrain += SIMD::size * 6;
        vonMisesStress  += SIMD::size;
        vonMisesStrain  += SIMD::size;

        for (int e = 0; e < element.elements; ++e) {
            ++enodes;
        }
    }

    template <typename Element>
    void element(Element &element)
    {
        constexpr double rgps = 1.0 / gps;
        SIMD eig[3], nu = element.ecf.poissonRatio[0];
        SIMD S[9] = {
                element.vS[0], element.vS[3], element.vS[5],
                element.vS[3], element.vS[1], element.vS[4],
                element.vS[5], element.vS[4], element.vS[2]
        };
        eigSym33Desc(S, eig);
        for (int e = 0; e < element.elements; ++e) {
            principalStress[3 * e + 0] += rgps * eig[0][e];
            principalStress[3 * e + 1] += rgps * eig[1][e];
            principalStress[3 * e + 2] += rgps * eig[2][e];

            componentStress[6 * e + 0] += rgps * S[0][e];
            componentStress[6 * e + 1] += rgps * S[4][e];
            componentStress[6 * e + 2] += rgps * S[8][e];
            componentStress[6 * e + 3] += rgps * S[1][e];
            componentStress[6 * e + 4] += rgps * S[5][e];
            componentStress[6 * e + 5] += rgps * S[2][e];

            vonMisesStress[e] += rgps * sqrt(.5 * ((eig[0][e] - eig[1][e]) * (eig[0][e] - eig[1][e]) + (eig[1][e] - eig[2][e]) * (eig[1][e] - eig[2][e]) + (eig[2][e] - eig[0][e]) * (eig[2][e] - eig[0][e])));
        }

        SIMD E[9] = {
                load1(.5) * element.C2[0] - load1(.5), load1(.5) * element.C2[1]            , load1(.5) * element.C2[2],
                load1(.5) * element.C2[3]            , load1(.5) * element.C2[4] - load1(.5), load1(.5) * element.C2[5],
                load1(.5) * element.C2[6]            , load1(.5) * element.C2[7]            , load1(.5) * element.C2[8] - load1(.5)
        };
        eigSym33Desc(E, eig);
        for (int e = 0; e < element.elements; ++e) {
            principalStrain[3 * e + 0] += rgps * eig[0][e];
            principalStrain[3 * e + 1] += rgps * eig[1][e];
            principalStrain[3 * e + 2] += rgps * eig[2][e];

            componentStrain[6 * e + 0] += rgps * E[0][e];
            componentStrain[6 * e + 1] += rgps * E[4][e];
            componentStrain[6 * e + 2] += rgps * E[8][e];
            componentStrain[6 * e + 3] += rgps * E[1][e];
            componentStrain[6 * e + 4] += rgps * E[5][e];
            componentStrain[6 * e + 5] += rgps * E[2][e];

            vonMisesStrain[e] += rgps * sqrt(.5 * ((eig[0][e] - eig[1][e]) * (eig[0][e] - eig[1][e]) + (eig[1][e] - eig[2][e]) * (eig[1][e] - eig[2][e]) + (eig[2][e] - eig[0][e]) * (eig[2][e] - eig[0][e]))) / (1 + nu[e]);
        }
    }

    template <typename Element>
    void node(Element &element, size_t n)
    {
        size_t index = 0;
        for (size_t nn = 0; nn < nodes; ++nn) {
            index += nn * element.NN[n][nn]; // get node index, only one values in NN is 1, other are 0
        }
        SIMD nu = element.ecf.poissonRatio[0];
        SIMD eigS[3], S[9] = {
                element.vS[0], element.vS[3], element.vS[5],
                element.vS[3], element.vS[1], element.vS[4],
                element.vS[5], element.vS[4], element.vS[2]
        };
        eigSym33Desc(S, eigS);
        SIMD eigE[3], E[9] = {
                load1(.5) * element.C2[0] - load1(.5), load1(.5) * element.C2[1]            , load1(.5) * element.C2[2],
                load1(.5) * element.C2[3]            , load1(.5) * element.C2[4] - load1(.5), load1(.5) * element.C2[5],
                load1(.5) * element.C2[6]            , load1(.5) * element.C2[7]            , load1(.5) * element.C2[8] - load1(.5)
        };
        eigSym33Desc(E, eigE);
        auto enodes = this->enodes;
        for (int e = 0; e < element.elements; ++e, ++enodes) {
            size_t node = enodes->at(index);
            double scale = multiplicity[node];
            principalStressAvg[3 * node + 0] += scale * eigS[0][e];
            principalStressAvg[3 * node + 1] += scale * eigS[1][e];
            principalStressAvg[3 * node + 2] += scale * eigS[2][e];
            componentStressAvg[6 * node + 0] += scale * S[0][e];
            componentStressAvg[6 * node + 1] += scale * S[4][e];
            componentStressAvg[6 * node + 2] += scale * S[8][e];
            componentStressAvg[6 * node + 3] += scale * S[1][e];
            componentStressAvg[6 * node + 4] += scale * S[5][e];
            componentStressAvg[6 * node + 5] += scale * S[2][e];
            vonMisesStressAvg [    node    ] += scale * sqrt(.5 * ((eigS[0][e] - eigS[1][e]) * (eigS[0][e] - eigS[1][e]) + (eigS[1][e] - eigS[2][e]) * (eigS[1][e] - eigS[2][e]) + (eigS[2][e] - eigS[0][e]) * (eigS[2][e] - eigS[0][e])));

            principalStrainAvg[3 * node + 0] += scale * eigE[0][e];
            principalStrainAvg[3 * node + 1] += scale * eigE[1][e];
            principalStrainAvg[3 * node + 2] += scale * eigE[2][e];
            componentStrainAvg[6 * node + 0] += scale * E[0][e];
            componentStrainAvg[6 * node + 1] += scale * E[4][e];
            componentStrainAvg[6 * node + 2] += scale * E[8][e];
            componentStrainAvg[6 * node + 3] += scale * E[1][e];
            componentStrainAvg[6 * node + 4] += scale * E[5][e];
            componentStrainAvg[6 * node + 5] += scale * E[2][e];
            vonMisesStrainAvg [    node    ] += scale * sqrt(.5 * ((eigE[0][e] - eigE[1][e]) * (eigE[0][e] - eigE[1][e]) + (eigE[1][e] - eigE[2][e]) * (eigE[1][e] - eigE[2][e]) + (eigE[2][e] - eigE[0][e]) * (eigE[2][e] - eigE[0][e]))) / (1 + nu[e]);
        }
    }

    void final()
    {
        SIMD X[13];

        SIMD eig[3], E[9] = {
                            X[ 6], load1(.5) * X[ 9], load1(.5) * X[11],
                load1(.5) * X[ 9],             X[ 7], load1(.5) * X[10],
                load1(.5) * X[11], load1(.5) * X[10],             X[ 8]
        };
        eigSym33Desc(E, eig);
//        SIMD intensity = max(abs(eig[0] - eig[1]), abs(eig[1] - eig[2]), abs(eig[0] - eig[2]));

        SIMD S[9] = {
                X[0], X[3], X[5],
                X[3], X[1], X[4],
                X[5], X[4], X[2]
        };
        eigSym33Desc(S, eig);
//        SIMD intensity = max(abs(eig[0] - eig[1]), abs(eig[1] - eig[2]), abs(eig[0] - eig[2]));

        SIMD vonMissesStress = sqrt(load1(.5) * ((eig[0] - eig[1]) * (eig[0] - eig[1]) + (eig[1] - eig[2]) * (eig[1] - eig[2]) + (eig[2] - eig[0]) * (eig[2] - eig[0])));
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_ */
