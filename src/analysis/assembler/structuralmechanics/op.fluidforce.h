
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FLUIDFORCE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FLUIDFORCE_H_

#include "analysis/assembler/general/subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "wrappers/simd/simd.h"

namespace espreso {

struct FluidForce: SubKernel {
    const char* name() const { return "FluidForces"; }

    serializededata<esint, esint>::const_iterator enodes, end;
    double * source;

    FluidForce()
    : enodes(info::mesh->elements->nodes->cbegin()),
      end(info::mesh->elements->nodes->cend()),
      source(nullptr)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION;
    }

    void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double * source)
    {
        this->enodes = enodes;
        this->end = end;
        this->source = source;
        this->isactive = 1;
    }
};

template <size_t ndim>
struct FluidForceKernel: FluidForce {
    FluidForceKernel(const FluidForce &base): FluidForce(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t d = 0; d < ndim; ++d) {
                element.f[d][s] = element.f[d][s] + source[ndim * enodes->at(0) + d];
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FLUIDFORCE_H_ */
