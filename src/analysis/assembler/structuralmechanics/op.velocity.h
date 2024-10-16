
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_VELOCITY_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_VELOCITY_H_

#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct Velocity: SubKernel {
    const char* name() const { return "Velocity"; }

    serializededata<esint, esint>::const_iterator enodes, end;
    double * target;

    Velocity()
    : enodes(info::mesh->elements->nodes->cbegin()),
      end(info::mesh->elements->nodes->cend()),
      target(nullptr)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
    }

    void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double * target)
    {
        this->enodes = enodes;
        this->end = end;
        this->target = target;
        this->isactive = 1;
    }
};

template <size_t nodes, size_t ndim>
struct InitialVelocityKernel: Velocity {
    InitialVelocityKernel(const Velocity &base): Velocity(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t n = 0; n < nodes; ++n) {
                for (size_t d = 0; d < ndim; ++d) {
                    target[enodes->at(n) * ndim + d] = element.velocity[n][d][s];
                }
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_VELOCITY_H_ */
