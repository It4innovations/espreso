
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

struct CoordinateRotation: public ActionOperator {
	CoordinateRotation(
			int interval,
			const ParameterData &coordinates,
			const ParameterData &rotation,
			ParameterData &angle)
	: coordinates(coordinates, interval),
	  rotation(rotation, interval),
	  angle(angle, interval)
	{

	}

	InputParameterIterator coordinates, rotation; // rotation -> angle for cartesian, center for cylinder and sphere
	OutputParameterIterator angle;

	void operator++()
	{
		++coordinates; ++rotation;
		++angle;
	}

	void move(int n)
	{
		coordinates += n;
		rotation += n;
		angle += n;
	}
};

template<size_t nodes, size_t gps>
struct CartesianRotation2D: CoordinateRotation {
	using CoordinateRotation::CoordinateRotation;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			angle[2 * gpindex + 0] = std::cos(M_PI * rotation[gpindex] / 180);
			angle[2 * gpindex + 1] = std::sin(M_PI * rotation[gpindex] / 180);
		}
	}
};

template<size_t nodes, size_t gps>
struct CartesianRotation3D: CoordinateRotation {
	using CoordinateRotation::CoordinateRotation;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			angle[3 * gpindex + 0] = std::cos(M_PI * rotation[3 * gpindex + 0] / 180);
			angle[3 * gpindex + 1] = std::cos(M_PI * rotation[3 * gpindex + 1] / 180);
			angle[3 * gpindex + 2] = std::cos(M_PI * rotation[3 * gpindex + 2] / 180);
			angle[3 * gpindex + 3] = std::sin(M_PI * rotation[3 * gpindex + 0] / 180);
			angle[3 * gpindex + 4] = std::sin(M_PI * rotation[3 * gpindex + 1] / 180);
			angle[3 * gpindex + 5] = std::sin(M_PI * rotation[3 * gpindex + 2] / 180);
		}
	}
};

template<size_t nodes, size_t gps>
struct CylindricalRotation2D: CoordinateRotation {
	using CoordinateRotation::CoordinateRotation;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double rot = std::atan2(coordinates[2 * gpindex + 1] - rotation[2 * gpindex + 1], coordinates[2 * gpindex + 0] - rotation[2 * gpindex + 0]);
			angle[2 * gpindex + 0] = std::cos(rot);
			angle[2 * gpindex + 1] = std::sin(rot);
		}
	}
};

template<size_t nodes, size_t gps>
struct CylindricalRotation3D: CoordinateRotation {
	using CoordinateRotation::CoordinateRotation;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double rot = std::atan2(coordinates[3 * gpindex + 1] - rotation[2 * gpindex + 1], coordinates[3 * gpindex + 0] - rotation[2 * gpindex + 0]);
			angle[3 * gpindex + 0] = 1;
			angle[3 * gpindex + 1] = 1;
			angle[3 * gpindex + 2] = std::cos(rot);
			angle[3 * gpindex + 3] = 0;
			angle[3 * gpindex + 4] = 0;
			angle[3 * gpindex + 5] = std::sin(rot);
		}
	}
};

template<size_t nodes, size_t gps>
struct SphericalRotation3D: CoordinateRotation {
	using CoordinateRotation::CoordinateRotation;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double x = coordinates[3 * gpindex + 0] - rotation[3 * gpindex + 0];
			double y = coordinates[3 * gpindex + 1] - rotation[3 * gpindex + 1];
			double z = coordinates[3 * gpindex + 2] - rotation[3 * gpindex + 2];
			double azimut = std::atan2(y, x);
			double r = std::sqrt(x * x + y * y + z * z);
			double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z * z + x * x), y);

			angle[3 * gpindex + 0] = 1;
			angle[3 * gpindex + 1] = std::cos(elevation);
			angle[3 * gpindex + 2] = std::cos(azimut);
			angle[3 * gpindex + 3] = 0;
			angle[3 * gpindex + 4] = std::sin(elevation);
			angle[3 * gpindex + 5] = std::sin(azimut);
		}
	}
};

struct ConductivityRotation: public ActionOperator {
	ConductivityRotation(
			int interval,
			const ParameterData &angle,
			ParameterData &result)
	: angle(angle, interval),
	  result(result, interval)
	{

	}

	InputParameterIterator angle;
	OutputParameterIterator result;

	void operator++()
	{
		++angle;
		++result;
	}

	void move(int n)
	{
		angle += n;
		result += n;
	}
};

template<size_t nodes, size_t gps>
struct ConductivityRotation2D: ConductivityRotation {
	using ConductivityRotation::ConductivityRotation;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			rotate2x2(angle[2 * gpindex + 0], angle[2 * gpindex + 1], result.data + 4 * gpindex);
		}
	}
};

template<size_t nodes, size_t gps>
struct ConductivityRotation3D: ConductivityRotation {
	using ConductivityRotation::ConductivityRotation;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			rotate3x3(angle.data + 6 * gpindex, angle.data + 6 * gpindex + 3, result.data + 9 * gpindex);
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_ */
