
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_CONDUCTIVITY_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_CONDUCTIVITY_H_

#include "analysis/assembler/general/op.coordinatesystem.h"
#include "config/ecf/material/thermalconductivity.h"
#include "config/ecf/material/coordinatesystem.h"

namespace espreso {

struct Conductivity: SubKernel {
    const char* name() const { return "ConductivityKernel"; }

    Conductivity()
    : conductivity(nullptr), coordinateSystem(nullptr), type(CoordinateSystemConfiguration::TYPE::CARTESIAN), rotated(false), constRotation(true)
    {
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(const ThermalConductivityConfiguration *conductivity, const CoordinateSystemConfiguration *coordinateSystem)
    {
        this->conductivity = conductivity;
        this->coordinateSystem = coordinateSystem;
        this->isconst = !(conductivity->needCoordinates() || conductivity->needTemperature()) && coordinateSystem->isConst();
        this->type = coordinateSystem->type;
        this->rotated = coordinateSystem->isRotated() && conductivity->model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
        this->constRotation = coordinateSystem->type == CoordinateSystemConfiguration::TYPE::CARTESIAN && coordinateSystem->isConst();
        this->isactive = !this->isconst || !this->constRotation;
    }

    const ThermalConductivityConfiguration *conductivity;
    const CoordinateSystemConfiguration *coordinateSystem;
    CoordinateSystemConfiguration::TYPE type;
    bool rotated, constRotation;
};

template <size_t ndim> struct ConductivityKernel;

template <> struct ConductivityKernel<2>: Conductivity {
    ConductivityKernel(const Conductivity &base): Conductivity(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        CoordinateSystem<2, 1>::simd(element, type, false);
        simd(element, 0);
    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        // 0 1
        // 2 3
        if (rotated) {
            switch (conductivity->model) {
            case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
            {

            } break;
            case ThermalConductivityConfiguration::MODEL::DIAGONAL:
            {
                CoordinateSystem<2, 1>::simd(element, type, constRotation);
                SIMD c00 = element.ecf.conductivity[0];
                SIMD c11 = element.ecf.conductivity[3];
                SIMD cos = element.rotation.cossin[0];
                SIMD sin = element.rotation.cossin[1];

                element.conductivity[0] = (cos * c00) * cos + (sin * c11) * sin;
                element.conductivity[1] = (cos * c00) * sin - (sin * c11) * cos;
                element.conductivity[3] = (sin * c00) * sin + (cos * c11) * cos;
            } break;
            case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
            {
                CoordinateSystem<2, 1>::simd(element, type, constRotation);
                SIMD c00 = element.ecf.conductivity[0];
                SIMD c01 = element.ecf.conductivity[1], c11 = element.ecf.conductivity[3];
                SIMD cos = element.rotation.cossin[0];
                SIMD sin = element.rotation.cossin[1];

                element.conductivity[0] = (cos * c00 - sin * c01) * cos - (cos * c01 - sin * c11) * sin;
                element.conductivity[1] = (cos * c00 - sin * c01) * sin + (cos * c01 - sin * c11) * cos;
                element.conductivity[3] = (sin * c00 + cos * c01) * sin + (sin * c01 + cos * c11) * cos;
            } break;
            case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
            {
                CoordinateSystem<2, 1>::simd(element, type, constRotation);
                SIMD c00 = element.ecf.conductivity[0], c10 = element.ecf.conductivity[2];
                SIMD c01 = element.ecf.conductivity[1], c11 = element.ecf.conductivity[3];
                SIMD cos = element.rotation.cossin[0];
                SIMD sin = element.rotation.cossin[1];

                element.conductivity[0] = (cos * c00 - sin * c10) * cos - (cos * c01 - sin * c11) * sin;
                element.conductivity[1] = (cos * c00 - sin * c10) * sin + (cos * c01 - sin * c11) * cos;
                element.conductivity[2] = (sin * c00 + cos * c10) * cos - (sin * c01 + cos * c11) * sin;
                element.conductivity[3] = (sin * c00 + cos * c10) * sin + (sin * c01 + cos * c11) * cos;
            } break;
            }
        } else {
            for (size_t i = 0; i < 4; ++i) {
                element.conductivity[i] = element.ecf.conductivity[i];
            }
        }
    }
};

template <> struct ConductivityKernel<3>: Conductivity {
    ConductivityKernel(const Conductivity &base): Conductivity(base) {}

    template <typename Element>
    void simd(Element &element)
    {
        CoordinateSystem<3, 1>::simd(element, type, false);
        simd(element, 0);
    }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        // 0 1 2
        // 3 4 5
        // 6 7 8
        if (rotated) {
            switch (conductivity->model) {
            case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
            {

            } break;
            case ThermalConductivityConfiguration::MODEL::DIAGONAL:
            {
                CoordinateSystem<3, 1>::simd(element, type, constRotation);
                SIMD cos0 = element.rotation.cossin[0];
                SIMD cos1 = element.rotation.cossin[1];
                SIMD cos2 = element.rotation.cossin[2];
                SIMD sin0 = element.rotation.cossin[3];
                SIMD sin1 = element.rotation.cossin[4];
                SIMD sin2 = element.rotation.cossin[5];

                SIMD t00 = cos1 * cos2;
                SIMD t01 = cos1 * sin2;
                SIMD t02 = -sin1;
                SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
                SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
                SIMD t12 = cos1 * sin0;
                SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
                SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
                SIMD t22 = cos0 * cos1;

                SIMD c00 = element.ecf.conductivity[0];
                SIMD c11 = element.ecf.conductivity[4];
                SIMD c22 = element.ecf.conductivity[8];

                SIMD a = t00 * c00, b = t10 * c11, c = t20 * c22;
                element.conductivity[0] = a * t00 + b * t10 + c * t20;
                element.conductivity[1] = a * t01 + b * t11 + c * t21;
                element.conductivity[2] = a * t02 + b * t12 + c * t22;

                a = t01 * c00, b = t11 * c11, c = t21 * c22;
                element.conductivity[4] = a * t01 + b * t11 + c * t21;
                element.conductivity[5] = a * t02 + b * t12 + c * t22;

                a = t02 * c00, b = t12 * c11, c = t22 * c22;
                element.conductivity[8] = a * t02 + b * t12 + c * t22;
            } break;
            case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
            {
                CoordinateSystem<3, 1>::simd(element, type, constRotation);
                SIMD cos0 = element.rotation.cossin[0];
                SIMD cos1 = element.rotation.cossin[1];
                SIMD cos2 = element.rotation.cossin[2];
                SIMD sin0 = element.rotation.cossin[3];
                SIMD sin1 = element.rotation.cossin[4];
                SIMD sin2 = element.rotation.cossin[5];

                SIMD t00 = cos1 * cos2;
                SIMD t01 = cos1 * sin2;
                SIMD t02 = -sin1;
                SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
                SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
                SIMD t12 = cos1 * sin0;
                SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
                SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
                SIMD t22 = cos0 * cos1;

                SIMD c00 = element.ecf.conductivity[0];
                SIMD c01 = element.ecf.conductivity[1], c11 = element.ecf.conductivity[4];
                SIMD c02 = element.ecf.conductivity[2], c12 = element.ecf.conductivity[5], c22 = element.ecf.conductivity[8];

                SIMD a = t00 * c00 + t10 * c01 + t20 * c02;
                SIMD b = t00 * c01 + t10 * c11 + t20 * c12;
                SIMD c = t00 * c02 + t10 * c12 + t20 * c22;
                element.conductivity[0] = a * t00 + b * t10 + c * t20;
                element.conductivity[1] = a * t01 + b * t11 + c * t21;
                element.conductivity[2] = a * t02 + b * t12 + c * t22;

                a = t01 * c00 + t11 * c01 + t21 * c02;
                b = t01 * c01 + t11 * c11 + t21 * c12;
                c = t01 * c02 + t11 * c12 + t21 * c22;
                element.conductivity[4] = a * t01 + b * t11 + c * t21;
                element.conductivity[5] = a * t02 + b * t12 + c * t22;

                a = t02 * c00 + t12 * c01 + t22 * c02;
                b = t02 * c01 + t12 * c11 + t22 * c12;
                c = t02 * c02 + t12 * c12 + t22 * c22;
                element.conductivity[8] = a * t02 + b * t12 + c * t22;
            } break;
            case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
            {
                CoordinateSystem<3, 1>::simd(element, type, constRotation);
                SIMD cos0 = element.rotation.cossin[0];
                SIMD cos1 = element.rotation.cossin[1];
                SIMD cos2 = element.rotation.cossin[2];
                SIMD sin0 = element.rotation.cossin[3];
                SIMD sin1 = element.rotation.cossin[4];
                SIMD sin2 = element.rotation.cossin[5];

                SIMD t00 = cos1 * cos2;
                SIMD t01 = cos1 * sin2;
                SIMD t02 = -sin1;
                SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
                SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
                SIMD t12 = cos1 * sin0;
                SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
                SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
                SIMD t22 = cos0 * cos1;

                SIMD c00 = element.ecf.conductivity[0], c10 = element.ecf.conductivity[3], c20 = element.ecf.conductivity[6];
                SIMD c01 = element.ecf.conductivity[1], c11 = element.ecf.conductivity[4], c21 = element.ecf.conductivity[7];
                SIMD c02 = element.ecf.conductivity[2], c12 = element.ecf.conductivity[5], c22 = element.ecf.conductivity[8];

                SIMD a = t00 * c00 + t10 * c10 + t20 * c20;
                SIMD b = t00 * c01 + t10 * c11 + t20 * c21;
                SIMD c = t00 * c02 + t10 * c12 + t20 * c22;
                element.conductivity[0] = a * t00 + b * t10 + c * t20;
                element.conductivity[1] = a * t01 + b * t11 + c * t21;
                element.conductivity[2] = a * t02 + b * t12 + c * t22;

                a = t01 * c00 + t11 * c10 + t21 * c20;
                b = t01 * c01 + t11 * c11 + t21 * c21;
                c = t01 * c02 + t11 * c12 + t21 * c22;
                element.conductivity[3] = a * t00 + b * t10 + c * t20;
                element.conductivity[4] = a * t01 + b * t11 + c * t21;
                element.conductivity[5] = a * t02 + b * t12 + c * t22;

                a = t02 * c00 + t12 * c10 + t22 * c20;
                b = t02 * c01 + t12 * c11 + t22 * c21;
                c = t02 * c02 + t12 * c12 + t22 * c22;
                element.conductivity[6] = a * t00 + b * t10 + c * t20;
                element.conductivity[7] = a * t01 + b * t11 + c * t21;
                element.conductivity[8] = a * t02 + b * t12 + c * t22;
            } break;
            }
        } else {
            for (size_t i = 0; i < 9; ++i) {
                element.conductivity[i] = element.ecf.conductivity[i];
            }
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_OP_CONDUCTIVITY_H_ */
