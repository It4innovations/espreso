
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_TEMPERATURE_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_TEMPERATURE_H_

#include "element.h"
#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct Temperature: SubKernel {
    const char* name() const { return "Temperature"; }

    serializededata<esint, esint>::const_iterator enodes, end;
    double * source;
    bool toGPs;

    Temperature()
    : enodes(info::mesh->elements->nodes->cbegin()),
      end(info::mesh->elements->nodes->cend()),
      source(nullptr), toGPs(false)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
    }

    void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double * source, bool toGPs)
    {
        this->enodes = enodes;
        this->end = end;
        this->source = source;
        this->toGPs = toGPs;
        this->isactive = 1;
    }
};

template <size_t nodes>
struct TemperatureKernel: Temperature {
    TemperatureKernel(const Temperature &base): Temperature(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t n = 0; n < nodes; ++n) {
                element.temperature.node[n][s] = source[enodes->at(n)];
            }
        }
    }
};

template <size_t nodes>
struct InitialTemperatureKernel: Temperature {
    InitialTemperatureKernel(const Temperature &base): Temperature(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t n = 0; n < nodes; ++n) {
                source[enodes->at(n)] = element.temperature.initial[n][s];
            }
        }
    }
};

template <size_t nodes>
struct TemperatureToGPsKernel: Temperature {
    TemperatureToGPsKernel(const Temperature &base): Temperature(base)
    {
        this->isactive = toGPs;
    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        element.temperature.gp = zeros();
        for (size_t n = 0; n < nodes; ++n) {
            element.temperature.gp = element.temperature.gp + load1(element.N[gp][n]) * element.temperature.node[n];
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_TEMPERATURE_H_ */
