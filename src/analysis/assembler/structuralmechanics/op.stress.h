
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_

#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/math.h"
#include "analysis/assembler/general/op.integration.h"
#include "analysis/assembler/structuralmechanics/op.material.h"
#include "config/ecf/physics/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nameddata.h"

namespace espreso {

struct Stress: SubKernel {
    const char* name() const { return "Stress"; }

    Stress()
    : enodes(info::mesh->elements->nodes->cbegin()),
      end(info::mesh->elements->nodes->cend()),
      target(nullptr),
      multiplicity(nullptr),
      principalStress(nullptr),
      componentStress(nullptr),
      vonMisesStress(nullptr),
      vonMisesStressEnd(nullptr)
    {
        isconst = false;
        action = SubKernel::SOLUTION;
    }

    serializededata<esint, esint>::const_iterator enodes, end;
    double *target, *multiplicity;
    double* principalStress, *componentStress, *vonMisesStress, *vonMisesStressEnd;

    void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double *target, double *multiplicity, size_t interval, NamedData *principalStress, NamedData *componentStress, NamedData *vonMisesStress)
    {
        this->enodes = enodes;
        this->end = end;
        this->target = target;
        this->multiplicity = multiplicity;
        this->principalStress = principalStress->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
        this->componentStress = componentStress->data.data() + 2 * info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
        this->vonMisesStress = vonMisesStress->data.data() + info::mesh->elements->eintervals[interval].begin;
        this->vonMisesStressEnd = vonMisesStress->data.data() + vonMisesStress->data.size();
        isactive = 1;
    }
};

template <size_t nodes, size_t ndim> struct StressKernel;

template <size_t nodes>
struct StressKernel<nodes, 2>: Stress {
    StressKernel(const Stress &base): Stress(base) {}

    template <typename Element>
    void simd(Element &element, size_t n)
    {
        // TODO
    }
};

template <size_t nodes>
struct StressKernel<nodes, 3>: Stress {
    StressKernel(const Stress &base): Stress(base) {}

    template <typename Element>
    void simd(Element &element, size_t n)
    {
        element.stress[ 0 * nodes + n] = element.vS[0];
        element.stress[ 1 * nodes + n] = element.vS[1];
        element.stress[ 2 * nodes + n] = element.vS[2];
        element.stress[ 3 * nodes + n] = element.vS[3];
        element.stress[ 4 * nodes + n] = element.vS[4];
        element.stress[ 5 * nodes + n] = element.vS[5];
        element.stress[ 6 * nodes + n] = load1(.5) * element.C2[0] - load1(.5);
        element.stress[ 7 * nodes + n] = load1(.5) * element.C2[4] - load1(.5);
        element.stress[ 8 * nodes + n] = load1(.5) * element.C2[8] - load1(.5);
        element.stress[ 9 * nodes + n] = load1(.5) * element.C2[1];
        element.stress[10 * nodes + n] = load1(.5) * element.C2[5];
        element.stress[11 * nodes + n] = load1(.5) * element.C2[2];
        element.stress[12 * nodes + n] = zeros();

//        printf("stress:");
//        for (int i = 0; i < 12; ++i) {
//            printf(" %+.5e", element.stress[i * nodes + n][0]);
//        }
//        printf("\n");
    }

    template <typename Element>
    void store(Element &element)
    {
        // TODO: compute normal in nodes
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t n = 0; n < nodes; ++n) {
                for (size_t d = 0; d < 3; ++d) {
                    target[3 * enodes->at(n) + d] += element.stress[d][s] * multiplicity[enodes->at(n)];
                }
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_STRESS_H_ */
