
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATESYSTEM_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATESYSTEM_H_

#include "config/ecf/material/coordinatesystem.h"

#include <cmath>

namespace espreso {

template <size_t ndim, size_t multiplicity> struct CoordinateSystemCartesian;
template <size_t ndim, size_t multiplicity> struct CoordinateSystemCylindric;

template <size_t multiplicity>
struct CoordinateSystemCartesian<2, multiplicity> {

    constexpr static double straightAngleRec = 1.0 / 180;

    template <typename Element>
    static void simd(Element &element)
    {
        SIMD angle = element.rotation.center[0];
        for (size_t s = 0; s < SIMD::size; ++s) {
            for (size_t m = 1; m <= multiplicity; ++m) {
                element.rotation.cossin[(m - 1) * 2 + 0][s] = std::cos(m * M_PI * angle[s] * straightAngleRec);
                element.rotation.cossin[(m - 1) * 2 + 1][s] = std::sin(m * M_PI * angle[s] * straightAngleRec);
            }
        }
    }
};

template <size_t multiplicity>
struct CoordinateSystemCartesian<3, multiplicity> {

    constexpr static double straightAngleRec = 1.0 / 180;

    template <typename Element>
    static void simd(Element &element)
    {
        SIMD angleX = element.rotation.center[0];
        SIMD angleY = element.rotation.center[1];
        SIMD angleZ = element.rotation.center[2];
        for (size_t s = 0; s < SIMD::size; ++s) {
            for (size_t m = 1; m <= multiplicity; ++m) {
                element.rotation.cossin[(m - 1) * 6 + 0][s] = std::cos(m * M_PI * angleX[s] * straightAngleRec);
                element.rotation.cossin[(m - 1) * 6 + 1][s] = std::cos(m * M_PI * angleY[s] * straightAngleRec);
                element.rotation.cossin[(m - 1) * 6 + 2][s] = std::cos(m * M_PI * angleZ[s] * straightAngleRec);
                element.rotation.cossin[(m - 1) * 6 + 3][s] = std::sin(m * M_PI * angleX[s] * straightAngleRec);
                element.rotation.cossin[(m - 1) * 6 + 4][s] = std::sin(m * M_PI * angleY[s] * straightAngleRec);
                element.rotation.cossin[(m - 1) * 6 + 5][s] = std::sin(m * M_PI * angleZ[s] * straightAngleRec);
            }
        }
    }
};

template <size_t multiplicity>
struct CoordinateSystemCylindric<2, multiplicity> {

    template <typename Element>
    static void simd(Element &element)
    {
        SIMD cooX =    element.coords.gp[0];
        SIMD cooY =    element.coords.gp[1];
        SIMD centerX = element.rotation.center[0];
        SIMD centerY = element.rotation.center[1];
        SIMD distanceX = cooX - centerX;
        SIMD distanceY = cooY - centerY;
        for (size_t s = 0; s < SIMD::size; ++s) {
            double rot = std::atan2(distanceY[s], distanceX[s]);
            for (size_t m = 1; m <= multiplicity; ++m) {
                element.rotation.cossin[(m - 1) * 2 + 0][s] = std::cos(m * rot);
                element.rotation.cossin[(m - 1) * 2 + 1][s] = std::sin(m * rot);
            }
        }
    }
};

template <size_t multiplicity>
struct CoordinateSystemCylindric<3, multiplicity> {

    template <typename Element>
    static void simd(Element &element)
    {
        SIMD cooX =    element.coords.gp[0];
        SIMD cooY =    element.coords.gp[1];
        SIMD centerX = element.rotation.center[0];
        SIMD centerY = element.rotation.center[1];
        SIMD distanceX = cooX - centerX;
        SIMD distanceY = cooY - centerY;
        for (size_t s = 0; s < SIMD::size; ++s) {
            double rot = std::atan2(distanceY[s], distanceX[s]);
            for (size_t m = 1; m <= multiplicity; ++m) {
                element.rotation.cossin[(m - 1) * 6 + 0][s] = 1;
                element.rotation.cossin[(m - 1) * 6 + 1][s] = 1;
                element.rotation.cossin[(m - 1) * 6 + 2][s] = std::cos(m * rot);
                element.rotation.cossin[(m - 1) * 6 + 3][s] = 0;
                element.rotation.cossin[(m - 1) * 6 + 4][s] = 0;
                element.rotation.cossin[(m - 1) * 6 + 5][s] = std::sin(m * rot);
            }
        }
    }
};

template <size_t multiplicity>
struct CoordinateSystemSpherical {

    template <typename Element>
    static void simd(Element &element)
    {
        SIMD x = element.coords.gp[0] - element.rotation.center[0];
        SIMD y = element.coords.gp[1] - element.rotation.center[1];
        SIMD z = element.coords.gp[2] - element.rotation.center[2];
        for (size_t s = 0; s < SIMD::size; ++s) {
            double azimut = std::atan2(y[s], x[s]);
            double r = std::sqrt(x[s] * x[s] + y[s] * y[s] + z[s] * z[s]);
            double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z[s] * z[s] + x[s] * x[s]), y[s]);
            for (size_t m = 1; m <= multiplicity; ++m) {
                element.rotation.cossin[(m - 1) * 6 + 0][s] = 1;
                element.rotation.cossin[(m - 1) * 6 + 1][s] = std::cos(m * elevation);
                element.rotation.cossin[(m - 1) * 6 + 2][s] = std::cos(m * azimut);
                element.rotation.cossin[(m - 1) * 6 + 3][s] = 0;
                element.rotation.cossin[(m - 1) * 6 + 4][s] = std::sin(m * elevation);
                element.rotation.cossin[(m - 1) * 6 + 5][s] = std::sin(m * azimut);
            }
        }
    }
};

template <size_t ndim, size_t multiplicity> struct CoordinateSystem;

template <size_t multiplicity> struct CoordinateSystem<2, multiplicity> {

    template <typename Element>
    static void simd(Element &element, CoordinateSystemConfiguration::TYPE type, int isconst)
    {
        switch (type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:   if (!isconst) CoordinateSystemCartesian<2, multiplicity>::simd(element); break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:               CoordinateSystemCylindric<2, multiplicity>::simd(element); break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:                                                                            break;
        }
    }
};

template <size_t multiplicity> struct CoordinateSystem<3, multiplicity> {

    template <typename Element>
    static void simd(Element &element, CoordinateSystemConfiguration::TYPE type, int isconst)
    {
        switch (type) {
        case CoordinateSystemConfiguration::TYPE::CARTESIAN:   if (!isconst) CoordinateSystemCartesian<3, multiplicity>::simd(element); break;
        case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:               CoordinateSystemCylindric<3, multiplicity>::simd(element); break;
        case CoordinateSystemConfiguration::TYPE::SPHERICAL:                 CoordinateSystemSpherical<   multiplicity>::simd(element); break;
        }
    }
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATESYSTEM_H_ */
