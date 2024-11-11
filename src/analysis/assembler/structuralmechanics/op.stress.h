
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_

#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/math.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nameddata.h"

namespace espreso {

struct Stress: SubKernel {
    const char* name() const { return "Stress"; }

    Stress()
    : principalStress(nullptr), componentStress(nullptr), vonMisesStress(nullptr), vonMisesStressEnd(nullptr)
    {
        isconst = false;
        action = SubKernel::SOLUTION;
    }

    double* principalStress, *componentStress, *vonMisesStress, *vonMisesStressEnd;

    void activate(size_t interval, NamedData *principalStress, NamedData *componentStress, NamedData *vonMisesStress)
    {
        this->principalStress = principalStress->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
        this->componentStress = componentStress->data.data() + 2 * info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
        this->vonMisesStress = vonMisesStress->data.data() + info::mesh->elements->eintervals[interval].begin;
        this->vonMisesStressEnd = vonMisesStress->data.data() + vonMisesStress->data.size();
        isactive = 1;
    }
};

template <size_t nodes, size_t gps, size_t ndim> struct StressKernel;

template <size_t nodes, size_t gps>
struct StressKernel<nodes, gps, 2>: Stress {
    StressKernel(const Stress &base): Stress(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        // TODO
    }
};

template <size_t nodes, size_t gps>
struct StressKernel<nodes, gps, 3>: Stress {
    StressKernel(const Stress &base): Stress(base) {}

    const double rgps = 1.0 / gps;

    template <typename Element>
    void simd(Element &element)
    {
        size_t size = std::min((size_t)SIMD::size, (size_t)(vonMisesStressEnd - vonMisesStress));

        SIMD scale = load1(rgps);
        SIMD CuB[9] = {
                element.sigma[0] * scale, element.sigma[3] * scale, element.sigma[5] * scale,
                element.sigma[3] * scale, element.sigma[1] * scale, element.sigma[4] * scale,
                element.sigma[5] * scale, element.sigma[4] * scale, element.sigma[2] * scale
        };

        double * __restrict__ component = componentStress;
        for (size_t s = 0; s < size; ++s) {
            component[6 * s + 0] = CuB[0][s];
            component[6 * s + 1] = CuB[4][s];
            component[6 * s + 2] = CuB[8][s];
            component[6 * s + 3] = CuB[1][s];
            component[6 * s + 4] = CuB[5][s];
            component[6 * s + 5] = CuB[2][s];
        }

        SIMD e[3];
        eigSym33Desc(CuB, e);

        double * __restrict__ principal = principalStress;
        for (size_t s = 0; s < size; ++s) {
            principal[3 * s + 0] = e[0][s];
            principal[3 * s + 1] = e[1][s];
            principal[3 * s + 2] = e[2][s];
        }

        SIMD vm = (e[0] - e[1]) * (e[0] - e[1]) + (e[0] - e[2]) * (e[0] - e[2]) + (e[1] - e[2]) * (e[1] - e[2]);
        vm = sqrt(load1(1. / 2) * vm);
        double * __restrict__ vonMises = vonMisesStress;
        for (size_t s = 0; s < size; ++s) {
            vonMises[s] = vm[s];
        }

        principalStress += 3 * SIMD::size;
        componentStress += 6 * SIMD::size;
        vonMisesStress  += 1 * SIMD::size;
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_ */
