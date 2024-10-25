
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FLUIDSTRESS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FLUIDSTRESS_H_

#include "analysis/assembler/general/subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "wrappers/simd/simd.h"

namespace espreso {

struct FluidStress: SubKernel {
    const char* name() const { return "FluidStress"; }

    serializededata<esint, esint>::const_iterator enodes, end;
    double * source;

    FluidStress()
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

template <size_t nodes, size_t ndim>
struct FluidStressGatherKernel: FluidStress {
    FluidStressGatherKernel(const FluidStress &base): FluidStress(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t n = 0; n < nodes; ++n) {
                for (size_t d = 0; d < ndim; ++d) {
                    element.coupling.stress[n][d][s] = source[ndim * enodes->at(n) + d];
                }
            }
        }
    }
};

template <size_t nodes, size_t ndim>
struct FluidStressKernel: FluidStress {
    FluidStressKernel(const FluidStress &base): FluidStress(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        for (size_t n = 0; n < nodes; ++n) {
            SIMD scale = element.det * load1(element.w[gp]) * load1(element.N[gp][n]);
            for (size_t d = 0; d < ndim; ++d) {
                element.f[d * nodes + n] = element.f[d * nodes + n] + scale * element.coupling.stress[n][d];
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FLUIDSTRESS_H_ */
