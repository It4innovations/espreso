
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FORCE_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FORCE_H_

#include "analysis/assembler/general/boundarycondition.h"

namespace espreso {

struct Force: public BoundaryCondition {
    const char* name() const { return "Force"; }

    Force()
    {

    }

    void activate(ECFExpressionVector *force)
    {
        BoundaryCondition::activate(force);
    }

    void activate(double *source, esint offset, esint size, serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end)
    {
        BoundaryCondition::activate(source, offset, size, enodes, end);
    }

    void activate(double *source, esint offset, esint size)
    {
        BoundaryCondition::activate(source, offset, size);
    }
};

template <size_t nodes, size_t ndim> struct ForceSetter: Force {
    const double scale = 1. / nodes;

    ForceSetter(const Force &base): Force(base)
    {
        isactive = this->source != nullptr;
    }

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t s = 0; s < SIMD::size && offset < size; ++s, ++offset) {
            for (size_t d = 0; d < ndim; ++d) {
                element.ecf.force[d][s] = scale * source[d + ndim * offset];
            }
        }
    }
};

template <size_t nodes, size_t ndim> struct ForceKernel: public Force {
    ForceKernel(const Force &base): Force(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t d = 0; d < ndim; ++d) {
            for (size_t n = 0; n < nodes; ++n) {
                element.f[d * nodes + n] = element.f[d * nodes + n] + element.ecf.force[d];
            }
        }
    }
};

template <size_t ndim>
struct ForceSetter<1, ndim>: Force {
    ForceSetter(const Force &base): Force(base)
    {
        isactive = this->source != nullptr;
    }

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t d = 0; d < ndim; ++d) {
                element.ecf.force[d][s] = source[ndim * enodes->at(0) + d];
            }
        }
    }
};



template <size_t ndim> struct ForceKernel<1, ndim>: public Force {
    ForceKernel(const Force &base): Force(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t d = 0; d < ndim; ++d) {
            element.f[d] = element.f[d] + element.ecf.force[d];
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OP_FORCE_H_ */
