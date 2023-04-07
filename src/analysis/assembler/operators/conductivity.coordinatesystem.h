
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_CONDUCTIVITY_COORDINATESYSTEM_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_CONDUCTIVITY_COORDINATESYSTEM_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/module/heattransfer.element.h"
#include "math/simd/simd.h"

#include <cmath>

namespace espreso {

struct HeatTransferCoordinateSystem: ActionOperator {
	const char* name() const { return "HeatTransferCoordinateSystem"; }
	HeatTransferCoordinateSystem(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}
};


template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferCoordinateSystemCartesian;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferCoordinateSystemCylindric;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferCoordinateSystemSpherical;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferCoordinateSystemApply;

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct HeatTransferCoordinateSystemCartesian<nodes, gps, 2, edim, etype, Physics>: HeatTransferCoordinateSystem, Physics {
	using HeatTransferCoordinateSystem::HeatTransferCoordinateSystem;

	constexpr static double straightAngleRec = 1.0 / 180;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angle = element.ecf.center[gp][0];
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.cossin[gp][0][s] = std::cos(M_PI * angle[s] * straightAngleRec);
				element.cossin[gp][1][s] = std::sin(M_PI * angle[s] * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct HeatTransferCoordinateSystemCartesian<nodes, gps, 3, edim, etype, Physics>: HeatTransferCoordinateSystem, Physics {
	using HeatTransferCoordinateSystem::HeatTransferCoordinateSystem;

	constexpr static double straightAngleRec = 1.0 / 180;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angleX = element.ecf.center[gp][0];
			SIMD angleY = element.ecf.center[gp][1];
			SIMD angleZ = element.ecf.center[gp][2];
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.cossin[gp][0][s] = std::cos(M_PI * angleX[s] * straightAngleRec);
				element.cossin[gp][1][s] = std::cos(M_PI * angleY[s] * straightAngleRec);
				element.cossin[gp][2][s] = std::cos(M_PI * angleZ[s] * straightAngleRec);
				element.cossin[gp][3][s] = std::sin(M_PI * angleX[s] * straightAngleRec);
				element.cossin[gp][4][s] = std::sin(M_PI * angleY[s] * straightAngleRec);
				element.cossin[gp][5][s] = std::sin(M_PI * angleZ[s] * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct HeatTransferCoordinateSystemCylindric<nodes, gps, 2, edim, etype, Physics>: HeatTransferCoordinateSystem, Physics {

	HeatTransferCoordinateSystemCylindric(size_t interval): HeatTransferCoordinateSystem(interval)
	{
		isconst = false;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cooX =    element.gpcoords[gp][0];
			SIMD cooY =    element.gpcoords[gp][1];
			SIMD centerX = element.ecf.center[gp][0];
			SIMD centerY = element.ecf.center[gp][1];
			SIMD distanceX = cooX - centerX;
			SIMD distanceY = cooY - centerY;
			for (size_t s = 0; s < SIMD::size; ++s) {
				double rot = std::atan2(distanceY[s], distanceX[s]);
				element.cossin[gp][0][s] = std::cos(rot);
				element.cossin[gp][1][s] = std::sin(rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct HeatTransferCoordinateSystemCylindric<nodes, gps, 3, edim, etype, Physics>: HeatTransferCoordinateSystem, Physics {

	HeatTransferCoordinateSystemCylindric(size_t interval): HeatTransferCoordinateSystem(interval)
	{
		isconst = false;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cooX =    element.gpcoords[gp][0];
			SIMD cooY =    element.gpcoords[gp][1];
			SIMD centerX = element.ecf.center[gp][0];
			SIMD centerY = element.ecf.center[gp][1];
			SIMD distanceX = cooX - centerX;
			SIMD distanceY = cooY - centerY;
			for (size_t s = 0; s < SIMD::size; ++s) {
				double rot = std::atan2(distanceY[s], distanceX[s]);
				element.cossin[gp][0][s] = 1;
				element.cossin[gp][1][s] = 1;
				element.cossin[gp][2][s] = std::cos(rot);
				element.cossin[gp][3][s] = 0;
				element.cossin[gp][4][s] = 0;
				element.cossin[gp][5][s] = std::sin(rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct HeatTransferCoordinateSystemSpherical<nodes, gps, 3, edim, etype, Physics>: HeatTransferCoordinateSystem, Physics {

	HeatTransferCoordinateSystemSpherical(size_t interval): HeatTransferCoordinateSystem(interval)
	{
		isconst = false;
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD x = element.gpcoords[gp][0] - element.ecf.center[gp][0];
			SIMD y = element.gpcoords[gp][1] - element.ecf.center[gp][1];
			SIMD z = element.gpcoords[gp][2] - element.ecf.center[gp][2];
			for (size_t s = 0; s < SIMD::size; ++s) {
				double azimut = std::atan2(y[s], x[s]);
				double r = std::sqrt(x[s] * x[s] + y[s] * y[s] + z[s] * z[s]);
				double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z[s] * z[s] + x[s] * x[s]), y[s]);
				element.cossin[gp][0][s] = 1;
				element.cossin[gp][1][s] = std::cos(elevation);
				element.cossin[gp][2][s] = std::cos(azimut);
				element.cossin[gp][3][s] = 0;
				element.cossin[gp][4][s] = std::sin(elevation);
				element.cossin[gp][5][s] = std::sin(azimut);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferCoordinateSystemApply<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferCoordinateSystem, Physics {
	using HeatTransferCoordinateSystem::HeatTransferCoordinateSystem;

	// |cos , -sin| |c[0] , c[2]| | cos , sin|
	// |sin ,  cos| |c[3] , c[1]| |-sin , cos|
	// |cos * c[0] - sin * c[3] , cos * c[2] - sin * c[1]|  | cos , sin|
	// |sin * c[0] + cos * c[3] , sin * c[2] + cos * c[1]|  |-sin , cos|

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.ecf.conductivity[gp][0];
			SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][2];
			SIMD cos = element.cossin[gp][0];
			SIMD sin = element.cossin[gp][1];

			element.conductivity[gp][0] = (cos * c00 - sin * c01) * cos - (cos * c01 - sin * c11) * sin;
			element.conductivity[gp][1] = (cos * c00 - sin * c01) * sin + (cos * c01 - sin * c11) * cos;
			element.conductivity[gp][2] = (sin * c00 + cos * c01) * sin + (sin * c01 + cos * c11) * cos;
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferCoordinateSystemApply<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: HeatTransferCoordinateSystem, Physics {
	using HeatTransferCoordinateSystem::HeatTransferCoordinateSystem;

	// |cos , -sin| |c[0] , c[2]| | cos , sin|
	// |sin ,  cos| |c[3] , c[1]| |-sin , cos|
	// |cos * c[0] - sin * c[3] , cos * c[2] - sin * c[1]|  | cos , sin|
	// |sin * c[0] + cos * c[3] , sin * c[2] + cos * c[1]|  |-sin , cos|

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.ecf.conductivity[gp][0], c10 = element.ecf.conductivity[gp][2];
			SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][3];
			SIMD cos = element.cossin[gp][0];
			SIMD sin = element.cossin[gp][1];

			element.conductivity[gp][0] = (cos * c00 - sin * c10) * cos - (cos * c01 - sin * c11) * sin;
			element.conductivity[gp][1] = (cos * c00 - sin * c10) * sin + (cos * c01 - sin * c11) * cos;
			element.conductivity[gp][2] = (sin * c00 + cos * c10) * cos - (sin * c01 + cos * c11) * sin;
			element.conductivity[gp][3] = (sin * c00 + cos * c10) * sin + (sin * c01 + cos * c11) * cos;
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferCoordinateSystemApply<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferCoordinateSystem, Physics {
	using HeatTransferCoordinateSystem::HeatTransferCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cos0 = element.cossin[gp][0];
			SIMD cos1 = element.cossin[gp][1];
			SIMD cos2 = element.cossin[gp][2];
			SIMD sin0 = element.cossin[gp][3];
			SIMD sin1 = element.cossin[gp][4];
			SIMD sin2 = element.cossin[gp][5];

			SIMD t00 = cos1 * cos2;
			SIMD t01 = cos1 * sin2;
			SIMD t02 = -sin1;
			SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
			SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
			SIMD t12 = cos1 * sin0;
			SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
			SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
			SIMD t22 = cos0 * cos1;

			SIMD c00 = element.ecf.conductivity[gp][0];
			SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][3];
			SIMD c02 = element.ecf.conductivity[gp][2], c12 = element.ecf.conductivity[gp][4], c22 = element.ecf.conductivity[gp][5];

			SIMD a = t00 * c00 + t10 * c01 + t20 * c02;
			SIMD b = t00 * c01 + t10 * c11 + t20 * c12;
			SIMD c = t00 * c02 + t10 * c12 + t20 * c22;
			element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

			a = t01 * c00 + t11 * c01 + t21 * c02;
			b = t01 * c01 + t11 * c11 + t21 * c12;
			c = t01 * c02 + t11 * c12 + t21 * c22;
			element.conductivity[gp][3] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][4] = a * t02 + b * t12 + c * t22;

			a = t02 * c00 + t12 * c01 + t22 * c02;
			b = t02 * c01 + t12 * c11 + t22 * c12;
			c = t02 * c02 + t12 * c12 + t22 * c22;
			element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferCoordinateSystemApply<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, Physics>: HeatTransferCoordinateSystem, Physics {
	using HeatTransferCoordinateSystem::HeatTransferCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cos0 = element.cossin[gp][0];
			SIMD cos1 = element.cossin[gp][1];
			SIMD cos2 = element.cossin[gp][2];
			SIMD sin0 = element.cossin[gp][3];
			SIMD sin1 = element.cossin[gp][4];
			SIMD sin2 = element.cossin[gp][5];

			SIMD t00 = cos1 * cos2;
			SIMD t01 = cos1 * sin2;
			SIMD t02 = -sin1;
			SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
			SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
			SIMD t12 = cos1 * sin0;
			SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
			SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
			SIMD t22 = cos0 * cos1;

			SIMD c00 = element.ecf.conductivity[gp][0], c10 = element.ecf.conductivity[gp][3], c20 = element.ecf.conductivity[gp][6];
			SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][4], c21 = element.ecf.conductivity[gp][7];
			SIMD c02 = element.ecf.conductivity[gp][2], c12 = element.ecf.conductivity[gp][5], c22 = element.ecf.conductivity[gp][8];

			SIMD a = t00 * c00 + t10 * c10 + t20 * c20;
			SIMD b = t00 * c01 + t10 * c11 + t20 * c21;
			SIMD c = t00 * c02 + t10 * c12 + t20 * c22;
			element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

			a = t01 * c00 + t11 * c10 + t21 * c20;
			b = t01 * c01 + t11 * c11 + t21 * c21;
			c = t01 * c02 + t11 * c12 + t21 * c22;
			element.conductivity[gp][3] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][4] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;

			a = t02 * c00 + t12 * c10 + t22 * c20;
			b = t02 * c01 + t12 * c11 + t22 * c21;
			c = t02 * c02 + t12 * c12 + t22 * c22;
			element.conductivity[gp][6] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][7] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][8] = a * t02 + b * t12 + c * t22;
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_CONDUCTIVITY_COORDINATESYSTEM_H_ */
