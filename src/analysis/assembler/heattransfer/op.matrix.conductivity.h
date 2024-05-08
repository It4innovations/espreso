
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OP_MATRIX_CONDUCTIVITY_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OP_MATRIX_CONDUCTIVITY_H_

#include "analysis/assembler/general/subkernel.h"
#include "config/ecf/material/thermalconductivity.h"
#include "math/primitives/matrix_info.h"
#include "wrappers/simd/simd.h"


namespace espreso {

struct MatrixConductivity: public SubKernel {
    const char* name() const { return "MatrixConductivity"; }

    ThermalConductivityConfiguration::MODEL model;

    MatrixConductivity()
    : model(ThermalConductivityConfiguration::MODEL::ANISOTROPIC)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(ThermalConductivityConfiguration::MODEL model)
    {
        this->isactive = 1;
        this->model = model;
    }
};

template <size_t nodes, size_t ndim> struct MatrixConductivityKernel;

template <size_t nodes>
struct MatrixConductivityKernel<nodes, 2>: MatrixConductivity {
    MatrixConductivityKernel(const MatrixConductivity &base): MatrixConductivity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (model) {
        case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
        {
            SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]) * element.conductivity[0];
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                element.K[(n * nodes + n)] = element.K[(n * nodes + n)] + scale * (nx * nx + ny * ny);
                for (size_t m = n + 1; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];
                    SIMD k = scale * (nx * mx + ny * my);

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + k;
                    element.K[(m * nodes + n)] = element.K[(m * nodes + n)] + k;
                }
            }
        } break;
        case ThermalConductivityConfiguration::MODEL::DIAGONAL:
        {
            SIMD c00 = element.conductivity[0];
            SIMD c11 = element.conductivity[3];
            SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]);
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                SIMD a = c00 * nx;
                SIMD b = c11 * ny;
                element.K[(n * nodes + n)] = element.K[(n * nodes + n)] + scale * (a * nx + b * ny);
                for (size_t m = n + 1; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];
                    SIMD k = scale * (a * mx + b * my);

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + k;
                    element.K[(m * nodes + n)] = element.K[(m * nodes + n)] + k;
                }
            }
        } break;
        case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
        {
            SIMD c00 = element.conductivity[0];
            SIMD c10 = element.conductivity[1], c11 = element.conductivity[3];
            SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]);
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                SIMD a = nx * c00 + ny * c10;
                SIMD b = nx * c10 + ny * c11;
                element.K[(n * nodes + n)] = element.K[(n * nodes + n)] + scale * (a * nx + b * ny);
                for (size_t m = n + 1; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];
                    SIMD k = scale * (a * mx + b * my);

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + k;
                    element.K[(m * nodes + n)] = element.K[(m * nodes + n)] + k;
                }
            }
        } break;
        case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
        {
            SIMD c00 = element.conductivity[0], c01 = element.conductivity[2];
            SIMD c10 = element.conductivity[1], c11 = element.conductivity[3];
            SIMD scale = element.thickness.gp * element.det * load1(element.w[gp]);
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                SIMD a = nx * c00 + ny * c01;
                SIMD b = nx * c10 + ny * c11;
                for (size_t m = 0; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + scale * (a * mx + b * my);
                }
            }
        } break;
        }
    }
};

template <size_t nodes>
struct MatrixConductivityKernel<nodes, 3>: MatrixConductivity {
    MatrixConductivityKernel(const MatrixConductivity &base): MatrixConductivity(base) {}

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        switch (model) {
        case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
        {
            SIMD scale = element.det * load1(element.w[gp]) * element.conductivity[0];
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                SIMD nz = element.dND[n][2];
                element.K[(n * nodes + n)] = element.K[(n * nodes + n)] + scale * (nx * nx + ny * ny + nz * nz);
                for (size_t m = n + 1; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];
                    SIMD mz = element.dND[m][2];
                    SIMD k = scale * (nx * mx + ny * my + nz * mz);

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + k;
                    element.K[(m * nodes + n)] = element.K[(m * nodes + n)] + k;
                }
            }
        } break;
        case ThermalConductivityConfiguration::MODEL::DIAGONAL:
        {
            SIMD c00 = element.conductivity[0];
            SIMD c11 = element.conductivity[4];
            SIMD c22 = element.conductivity[8];
            SIMD scale = element.det * load1(element.w[gp]);
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                SIMD nz = element.dND[n][2];
                SIMD a = c00 * nx;
                SIMD b = c11 * ny;
                SIMD c = c22 * nz;
                element.K[(n * nodes + n)] = element.K[(n * nodes + n)] + scale * (a * nx + b * ny + c * nz);
                for (size_t m = n + 1; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];
                    SIMD mz = element.dND[m][2];
                    SIMD k = scale * (a * mx + b * my + c * mz);

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + k;
                    element.K[(m * nodes + n)] = element.K[(m * nodes + n)] + k;
                }
            }
        } break;
        case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
        {
            SIMD c00 = element.conductivity[0];
            SIMD c10 = element.conductivity[1], c11 = element.conductivity[4];
            SIMD c20 = element.conductivity[2], c21 = element.conductivity[5], c22 = element.conductivity[8];
            SIMD scale = element.det * load1(element.w[gp]);
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                SIMD nz = element.dND[n][2];
                SIMD a = nx * c00 + ny * c10 + nz * c20;
                SIMD b = nx * c10 + ny * c11 + nz * c21;
                SIMD c = nx * c20 + ny * c21 + nz * c22;
                element.K[(n * nodes + n)] = element.K[(n * nodes + n)] + scale * (a * nx + b * ny + c * nz);
                for (size_t m = n + 1; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];
                    SIMD mz = element.dND[m][2];
                    SIMD k = scale * (a * mx + b * my + c * mz);

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + k;
                    element.K[(m * nodes + n)] = element.K[(m * nodes + n)] + k;
                }
            }
        } break;
        case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
        {
            SIMD c00 = element.conductivity[0], c01 = element.conductivity[3], c02 = element.conductivity[6];
            SIMD c10 = element.conductivity[1], c11 = element.conductivity[4], c12 = element.conductivity[7];
            SIMD c20 = element.conductivity[2], c21 = element.conductivity[5], c22 = element.conductivity[8];
            SIMD scale = element.det * load1(element.w[gp]);
            for (size_t n = 0; n < nodes; ++n) {
                SIMD nx = element.dND[n][0];
                SIMD ny = element.dND[n][1];
                SIMD nz = element.dND[n][2];
                SIMD a = nx * c00 + ny * c01 + nz * c02;
                SIMD b = nx * c10 + ny * c11 + nz * c12;
                SIMD c = nx * c20 + ny * c21 + nz * c22;
                for (size_t m = 0; m < nodes; ++m) {
                    SIMD mx = element.dND[m][0];
                    SIMD my = element.dND[m][1];
                    SIMD mz = element.dND[m][2];

                    element.K[(n * nodes + m)] = element.K[(n * nodes + m)] + scale * (a * mx + b * my + c * mz);
                }
            }
        } break;
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OP_MATRIX_CONDUCTIVITY_H_ */
