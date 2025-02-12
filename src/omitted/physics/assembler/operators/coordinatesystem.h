
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/math.hpp"

namespace espreso {

struct CoordinateSystem: public Operator {
    CoordinateSystem(
            const ParameterData &coordinates,
            const ParameterData &rotation,
            ParameterData &matrix,
            int interval,
            int update = 1)
    : Operator(interval, matrix.isconst[interval], update && matrix.update[interval]),
      coordinates(coordinates, interval),
      rotation(rotation, interval),
      matrix(matrix, interval)
    {

    }

    InputParameterIterator coordinates, rotation;
    OutputParameterIterator matrix;

    void operator++()
    {
        ++coordinates; ++rotation;
        ++matrix;
    }
};

struct Cartesian2DCoordinateSystem: CoordinateSystem {
    Cartesian2DCoordinateSystem(
            const ParameterData &coordinates,
            const ElementExternalParameter<egps> &rotation,
            ParameterData &matrix,
            int interval): CoordinateSystem(coordinates, rotation, matrix, interval, rotation.isset[interval]), fixedRotation(rotation.isconst[interval])
    {
        angle();
    }

    double sin, cos;
    const bool fixedRotation;

    void angle()
    {
        sin = std::sin(M_PI * rotation[0] / 180);
        cos = std::cos(M_PI * rotation[0] / 180);
    }

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        if (!fixedRotation) {
            angle();
        }
        rotate2x2(cos, sin, matrix.data + 4 * gpindex);
    }
};

struct Cylindrical2DCoordinateSystem: CoordinateSystem {
    using CoordinateSystem::CoordinateSystem;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double angle = std::atan2(coordinates[2 * gpindex + 1] - rotation[2 * gpindex + 1], coordinates[2 * gpindex + 0] - rotation[2 * gpindex + 0]);
        double cos = std::cos(angle);
        double sin = std::sin(angle);
        rotate2x2(cos, sin, matrix.data + 4 * gpindex);
    }
};

struct Cartesian3DCoordinateSystem: CoordinateSystem {
    Cartesian3DCoordinateSystem(
            const ParameterData &coordinates,
            const ElementExternalParameter<ndim * egps> &rotation,
            ParameterData &matrix,
            int interval): CoordinateSystem(coordinates, rotation, matrix, interval, rotation.isset[interval]), fixedRotation(rotation.isconst[interval])
    {
        angle();
    }

    double sin[3], cos[3];
    const bool fixedRotation;

    void angle()
    {
        sin[0] = std::sin(M_PI * rotation[0] / 180);
        sin[1] = std::sin(M_PI * rotation[1] / 180);
        sin[2] = std::sin(M_PI * rotation[2] / 180);
        cos[0] = std::cos(M_PI * rotation[0] / 180);
        cos[1] = std::cos(M_PI * rotation[1] / 180);
        cos[2] = std::cos(M_PI * rotation[2] / 180);
    }

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        if (!fixedRotation) {
            angle();
        }
        rotate3x3(cos, sin, matrix.data + 9 * gpindex);
    }
};

struct Cylindrical3DCoordinateSystem: CoordinateSystem {
    using CoordinateSystem::CoordinateSystem;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double angle = std::atan2(coordinates[3 * gpindex + 1] - rotation[2 * gpindex + 1], coordinates[3 * gpindex + 0] - rotation[2 * gpindex + 0]);
        double cos[3] = { 1, 1, std::cos(angle) };
        double sin[3] = { 0, 0, std::sin(angle) };
        rotate3x3(cos, sin, matrix.data + 9 * gpindex);
    }
};

struct Spherical3DCoordinateSystem: CoordinateSystem {
    using CoordinateSystem::CoordinateSystem;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double x = coordinates[3 * gpindex + 0] - rotation[3 * gpindex + 0];
        double y = coordinates[3 * gpindex + 1] - rotation[3 * gpindex + 1];
        double z = coordinates[3 * gpindex + 2] - rotation[3 * gpindex + 2];
        double azimut = std::atan2(y, x);
        double r = std::sqrt(x * x + y * y + z * z);
        double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z * z + x * x), y);

        double cos[3] = { 1, std::cos(elevation), std::cos(azimut) };
        double sin[3] = { 0, std::sin(elevation), std::sin(azimut) };

        rotate3x3(cos, sin, matrix.data + 9 * gpindex);
    }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_ */
