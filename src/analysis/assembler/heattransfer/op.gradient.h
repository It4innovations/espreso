
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_GRADIENT_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_GRADIENT_H_

#include "analysis/assembler/general/subkernel.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nameddata.h"

namespace espreso {

struct TemperatureGradient: SubKernel {
    const char* name() const { return "TemperatureGradient"; }

    double *gradient, *end;

    TemperatureGradient()
    : gradient(nullptr), end(nullptr)
    {
        isconst = false;
        action = SubKernel::SOLUTION;
    }

    void activate(size_t interval, NamedData *gradient)
    {
        this->gradient = gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
        this->end = gradient->data.data() + gradient->data.size();
        isactive = 1;
    }
};

template <size_t nodes, size_t gps, size_t ndim> struct TemperatureGradientKernel;

template <size_t nodes, size_t gps>
struct TemperatureGradientKernel<nodes, gps, 2>: TemperatureGradient {
    TemperatureGradientKernel(const TemperatureGradient &base): TemperatureGradient(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        for (size_t n = 0; n < nodes; ++n) {
            element.gradient[0] = element.gradient[0] + element.dND[n][0] * element.temperature.node[n];
            element.gradient[1] = element.gradient[1] + element.dND[n][1] * element.temperature.node[n];
        }
    }

    template <typename Element>
    void store(Element &element)
    {
        SIMD scale = load1(1. / gps);
        element.gradient[0] = element.gradient[0] * scale;
        element.gradient[1] = element.gradient[1] * scale;
        size_t size = std::min((size_t)SIMD::size, (size_t)(end - gradient));
        double * __restrict__ out = gradient;
        for (size_t s = 0; s < size; ++s) {
            out[2 * s + 0] = element.gradient[0][s];
            out[2 * s + 1] = element.gradient[1][s];
        }
        element.gradient[0] = element.gradient[1] = zeros();
        gradient += 2 * size;
    }
};

template <size_t nodes, size_t gps>
struct TemperatureGradientKernel<nodes, gps, 3>: TemperatureGradient {
    TemperatureGradientKernel(const TemperatureGradient &base): TemperatureGradient(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        for (size_t n = 0; n < nodes; ++n) {
            element.gradient[0] = element.gradient[0] + element.dND[n][0] * element.temperature.node[n];
            element.gradient[1] = element.gradient[1] + element.dND[n][1] * element.temperature.node[n];
            element.gradient[2] = element.gradient[2] + element.dND[n][2] * element.temperature.node[n];
        }
    }

    template <typename Element>
    void store(Element &element)
    {
        SIMD scale = load1(1. / gps);
        element.gradient[0] = element.gradient[0] * scale;
        element.gradient[1] = element.gradient[1] * scale;
        element.gradient[2] = element.gradient[2] * scale;
        size_t size = std::min((size_t)SIMD::size, (size_t)(end - gradient));
        double * __restrict__ out = gradient;
        for (size_t s = 0; s < size; ++s) {
            out[3 * s + 0] = element.gradient[0][s];
            out[3 * s + 1] = element.gradient[1][s];
            out[3 * s + 2] = element.gradient[2][s];
        }
        element.gradient[0] = element.gradient[1] = element.gradient[2] = zeros();
        gradient += 3 * size;
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_GRADIENT_H_ */
