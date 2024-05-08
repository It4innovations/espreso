
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_FLUX_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_FLUX_H_

#include "analysis/assembler/general/subkernel.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nameddata.h"
#include "config/ecf/material/thermalconductivity.h"

namespace espreso {

struct TemperatureFlux: SubKernel {
    const char* name() const { return "TemperatureFlux"; }

    double *flux, *end;
    ThermalConductivityConfiguration::MODEL model;

    TemperatureFlux()
    : flux(nullptr), end(nullptr), model(ThermalConductivityConfiguration::MODEL::ANISOTROPIC)
    {
        isconst = false;
        action = SubKernel::SOLUTION;
    }

    void activate(size_t interval, NamedData *flux, ThermalConductivityConfiguration::MODEL model)
    {
        this->flux = flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
        this->end = flux->data.data() + flux->data.size();
        this->model = model;
        isactive = 1;
    }
};

template <size_t nodes, size_t gps, size_t ndim> struct TemperatureFluxKernel;

template <size_t nodes, size_t gps>
struct TemperatureFluxKernel<nodes, gps, 2>: TemperatureFlux {
    TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (model) {
        case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0];
            element.flux[1] = element.flux[1] + fgp1 * element.conductivity[0];
        } break;
        case ThermalConductivityConfiguration::MODEL::DIAGONAL:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0];
            element.flux[1] = element.flux[1] + fgp1 * element.conductivity[3];
        } break;
        case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0] + fgp1 * element.conductivity[1];
            element.flux[1] = element.flux[1] + fgp0 * element.conductivity[1] + fgp1 * element.conductivity[3];
        } break;
        case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0] + fgp1 * element.conductivity[1];
            element.flux[1] = element.flux[1] + fgp0 * element.conductivity[2] + fgp1 * element.conductivity[3];
        } break;
        }
    }

    template <typename Element>
    void store(Element &element)
    {
        SIMD scale = load1(1. / gps);
        element.flux[0] = element.flux[0] * scale;
        element.flux[1] = element.flux[1] * scale;
        size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
        double * __restrict__ out = flux;
        for (size_t s = 0; s < size; ++s) {
            out[2 * s + 0] = element.flux[0][s];
            out[2 * s + 1] = element.flux[1][s];
        }
        element.flux[0] = element.flux[1] = zeros();
        flux += 2 * size;
    }
};

template <size_t nodes, size_t gps>
struct TemperatureFluxKernel<nodes, gps, 3>: TemperatureFlux {
    TemperatureFluxKernel(const TemperatureFlux &base): TemperatureFlux(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (model) {
        case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
                fgp2 = fgp2 + element.dND[n][2] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0];
            element.flux[1] = element.flux[1] + fgp1 * element.conductivity[0];
            element.flux[2] = element.flux[2] + fgp2 * element.conductivity[0];
        } break;
        case ThermalConductivityConfiguration::MODEL::DIAGONAL:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
                fgp2 = fgp2 + element.dND[n][2] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0];
            element.flux[1] = element.flux[1] + fgp1 * element.conductivity[4];
            element.flux[2] = element.flux[2] + fgp2 * element.conductivity[8];
        } break;
        case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
                fgp2 = fgp2 + element.dND[n][2] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0] + fgp1 * element.conductivity[1] + fgp2 * element.conductivity[2];
            element.flux[1] = element.flux[1] + fgp0 * element.conductivity[1] + fgp1 * element.conductivity[4] + fgp2 * element.conductivity[5];
            element.flux[2] = element.flux[2] + fgp0 * element.conductivity[2] + fgp1 * element.conductivity[5] + fgp2 * element.conductivity[8];
        } break;
        case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
        {
            SIMD fgp0 = zeros(), fgp1 = zeros(), fgp2 = zeros();
            for (size_t n = 0; n < nodes; ++n) {
                fgp0 = fgp0 + element.dND[n][0] * element.temperature.node[n];
                fgp1 = fgp1 + element.dND[n][1] * element.temperature.node[n];
                fgp2 = fgp2 + element.dND[n][2] * element.temperature.node[n];
            }
            element.flux[0] = element.flux[0] + fgp0 * element.conductivity[0] + fgp1 * element.conductivity[1] + fgp2 * element.conductivity[2];
            element.flux[1] = element.flux[1] + fgp0 * element.conductivity[3] + fgp1 * element.conductivity[4] + fgp2 * element.conductivity[5];
            element.flux[2] = element.flux[2] + fgp0 * element.conductivity[6] + fgp1 * element.conductivity[7] + fgp2 * element.conductivity[8];
        } break;
        }
    }

    template <typename Element>
    void store(Element &element)
    {
        SIMD scale = load1(1. / gps);
        element.flux[0] = element.flux[0] * scale;
        element.flux[1] = element.flux[1] * scale;
        element.flux[2] = element.flux[2] * scale;
        size_t size = std::min((size_t)SIMD::size, (size_t)(end - flux));
        double * __restrict__ out = flux;
        for (size_t s = 0; s < size; ++s) {
            out[3 * s + 0] = element.flux[0][s];
            out[3 * s + 1] = element.flux[1][s];
            out[3 * s + 2] = element.flux[2][s];
        }
        element.flux[0] = element.flux[1] = element.flux[2] = zeros();
        flux += 3 * size;
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_FLUX_H_ */
