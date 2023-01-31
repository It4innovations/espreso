
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_

#include "analysis/assembler/operator.h"

#include <cmath>

namespace espreso {

struct CoordinateSystem: ActionOperator {
	CoordinateSystem(size_t interval) {}
};


template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct CoordinateSystemCartesian;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct CoordinateSystemCylindric;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct CoordinateSystemSpherical;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct CoordinateSystemApply;

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct CoordinateSystemCartesian<nodes, gps, 2, edim, etype, Physics>: CoordinateSystem, Physics {
	using CoordinateSystem::CoordinateSystem;

	constexpr static double straightAngleRec = 1.0 / 180;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double angle = element.ecf.angle[gp][0];
			element.ecf.angle[gp][0] = std::cos(M_PI * angle * straightAngleRec);
			element.ecf.angle[gp][1] = std::sin(M_PI * angle * straightAngleRec);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angle = element.ecf.angle[gp][0];
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.angle[gp][0][s] = std::cos(M_PI * angle[s] * straightAngleRec);
				element.ecf.angle[gp][1][s] = std::sin(M_PI * angle[s] * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct CoordinateSystemCartesian<nodes, gps, 3, edim, etype, Physics>: CoordinateSystem, Physics {
	using CoordinateSystem::CoordinateSystem;

	constexpr static double straightAngleRec = 1.0 / 180;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double angleX = element.ecf.angle[gp][0];
			double angleY = element.ecf.angle[gp][1];
			double angleZ = element.ecf.angle[gp][2];
			element.ecf.angle[gp][0] = std::cos(M_PI * angleX * straightAngleRec);
			element.ecf.angle[gp][1] = std::cos(M_PI * angleY * straightAngleRec);
			element.ecf.angle[gp][2] = std::cos(M_PI * angleZ * straightAngleRec);
			element.ecf.angle[gp][3] = std::sin(M_PI * angleX * straightAngleRec);
			element.ecf.angle[gp][4] = std::sin(M_PI * angleY * straightAngleRec);
			element.ecf.angle[gp][5] = std::sin(M_PI * angleZ * straightAngleRec);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angleX = element.ecf.angle[gp][0];
			SIMD angleY = element.ecf.angle[gp][1];
			SIMD angleZ = element.ecf.angle[gp][2];
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.angle[gp][0][s] = std::cos(M_PI * angleX[s] * straightAngleRec);
				element.ecf.angle[gp][1][s] = std::cos(M_PI * angleY[s] * straightAngleRec);
				element.ecf.angle[gp][2][s] = std::cos(M_PI * angleZ[s] * straightAngleRec);
				element.ecf.angle[gp][3][s] = std::sin(M_PI * angleX[s] * straightAngleRec);
				element.ecf.angle[gp][4][s] = std::sin(M_PI * angleY[s] * straightAngleRec);
				element.ecf.angle[gp][5][s] = std::sin(M_PI * angleZ[s] * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct CoordinateSystemCylindric<nodes, gps, 2, edim, etype, Physics>: CoordinateSystem, Physics {

	CoordinateSystemCylindric(size_t interval): CoordinateSystem(interval)
	{
		isconst = false;
	}

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double cooX =    element.gpcoords[gp][0];
			double cooY =    element.gpcoords[gp][1];
			double centerX = element.ecf.angle[gp][0];
			double centerY = element.ecf.angle[gp][1];
			double rot = std::atan2(cooY - centerY, cooX - centerX);
			element.ecf.angle[gp][0] = std::cos(rot);
			element.ecf.angle[gp][1] = std::sin(rot);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cooX =    element.gpcoords[gp][0];
			SIMD cooY =    element.gpcoords[gp][1];
			SIMD centerX = element.ecf.angle[gp][0];
			SIMD centerY = element.ecf.angle[gp][1];
			SIMD distanceX = cooX - centerX;
			SIMD distanceY = cooY - centerY;
			for (size_t s = 0; s < SIMD::size; ++s) {
				double rot = std::atan2(distanceY[s], distanceX[s]);
				element.ecf.angle[gp][0][s] = std::cos(rot);
				element.ecf.angle[gp][1][s] = std::sin(rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct CoordinateSystemCylindric<nodes, gps, 3, edim, etype, Physics>: CoordinateSystem, Physics {

	CoordinateSystemCylindric(size_t interval): CoordinateSystem(interval)
	{
		isconst = false;
	}

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double cooX =    element.gpcoords[gp][0];
			double cooY =    element.gpcoords[gp][1];
			double centerX = element.ecf.angle[gp][0];
			double centerY = element.ecf.angle[gp][1];
			double rot = std::atan2(cooY - centerY, cooX - centerX);
			element.ecf.angle[gp][0] = 1;
			element.ecf.angle[gp][1] = 1;
			element.ecf.angle[gp][2] = std::cos(rot);
			element.ecf.angle[gp][3] = 0;
			element.ecf.angle[gp][4] = 0;
			element.ecf.angle[gp][5] = std::sin(rot);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cooX =    element.gpcoords[gp][0];
			SIMD cooY =    element.gpcoords[gp][1];
			SIMD centerX = element.ecf.angle[gp][0];
			SIMD centerY = element.ecf.angle[gp][1];
			SIMD distanceX = cooX - centerX;
			SIMD distanceY = cooY - centerY;
			for (size_t s = 0; s < SIMD::size; ++s) {
				double rot = std::atan2(distanceY[s], distanceX[s]);
				element.ecf.angle[gp][0][s] = 1;
				element.ecf.angle[gp][1][s] = 1;
				element.ecf.angle[gp][2][s] = std::cos(rot);
				element.ecf.angle[gp][3][s] = 0;
				element.ecf.angle[gp][4][s] = 0;
				element.ecf.angle[gp][5][s] = std::sin(rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct CoordinateSystemSpherical<nodes, gps, 3, edim, etype, Physics>: CoordinateSystem, Physics {

	CoordinateSystemSpherical(size_t interval): CoordinateSystem(interval)
	{
		isconst = false;
	}

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double x = element.gpcoords[gp][0] - element.ecf.angle[gp][0];
			double y = element.gpcoords[gp][1] - element.ecf.angle[gp][1];
			double z = element.gpcoords[gp][2] - element.ecf.angle[gp][2];
			double azimut = std::atan2(y, x);
			double r = std::sqrt(x * x + y * y + z * z);
			double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z * z + x * x), y);
			element.ecf.angle[gp][0] = 1;
			element.ecf.angle[gp][1] = std::cos(elevation);
			element.ecf.angle[gp][2] = std::cos(azimut);
			element.ecf.angle[gp][3] = 0;
			element.ecf.angle[gp][4] = std::sin(elevation);
			element.ecf.angle[gp][5] = std::sin(azimut);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD x = element.gpcoords[gp][0] - element.ecf.angle[gp][0];
			SIMD y = element.gpcoords[gp][1] - element.ecf.angle[gp][1];
			SIMD z = element.gpcoords[gp][2] - element.ecf.angle[gp][2];
			for (size_t s = 0; s < SIMD::size; ++s) {
				double azimut = std::atan2(y[s], x[s]);
				double r = std::sqrt(x[s] * x[s] + y[s] * y[s] + z[s] * z[s]);
				double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z[s] * z[s] + x[s] * x[s]), y[s]);
				element.ecf.angle[gp][0][s] = 1;
				element.ecf.angle[gp][1][s] = std::cos(elevation);
				element.ecf.angle[gp][2][s] = std::cos(azimut);
				element.ecf.angle[gp][3][s] = 0;
				element.ecf.angle[gp][4][s] = std::sin(elevation);
				element.ecf.angle[gp][5][s] = std::sin(azimut);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct CoordinateSystemApply<nodes, gps, 2, edim, etype, Physics>: CoordinateSystem, Physics {
	using CoordinateSystem::CoordinateSystem;

	// |cos , -sin| |c[0] , c[2]| | cos , sin|
	// |sin ,  cos| |c[3] , c[1]| |-sin , cos|
	// |cos * c[0] - sin * c[3] , cos * c[2] - sin * c[1]|  | cos , sin|
	// |sin * c[0] + cos * c[3] , sin * c[2] + cos * c[1]|  |-sin , cos|
	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double origin[4] = {
					element.ecf.conductivity[gp][0], element.ecf.conductivity[gp][1],
					element.ecf.conductivity[gp][2], element.ecf.conductivity[gp][3] };
			double cos = element.ecf.angle[gp][0];
			double sin = element.ecf.angle[gp][1];
			element.conductivity[gp][0] = (cos * origin[0] - sin * origin[2]) * cos - (cos * origin[1] - sin * origin[3]) * sin;
			element.conductivity[gp][1] = (cos * origin[0] - sin * origin[2]) * sin + (cos * origin[1] - sin * origin[3]) * cos;
			element.conductivity[gp][2] = (sin * origin[0] + cos * origin[2]) * cos - (sin * origin[1] + cos * origin[3]) * sin;
			element.conductivity[gp][3] = (sin * origin[0] + cos * origin[2]) * sin + (sin * origin[1] + cos * origin[3]) * cos;
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD origin0 = element.ecf.conductivity[gp][0];
			SIMD origin1 = element.ecf.conductivity[gp][1];
			SIMD origin2 = element.ecf.conductivity[gp][2];
			SIMD origin3 = element.ecf.conductivity[gp][3];
			SIMD cos = element.ecf.angle[gp][0];
			SIMD sin = element.ecf.angle[gp][1];

			element.conductivity[gp][0] = (cos * origin0 - sin * origin2) * cos - (cos * origin1 - sin * origin3) * sin;
			element.conductivity[gp][1] = (cos * origin0 - sin * origin2) * sin + (cos * origin1 - sin * origin3) * cos;
			element.conductivity[gp][2] = (sin * origin0 + cos * origin2) * cos - (sin * origin1 + cos * origin3) * sin;
			element.conductivity[gp][3] = (sin * origin0 + cos * origin2) * sin + (sin * origin1 + cos * origin3) * cos;
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct CoordinateSystemApply<nodes, gps, 3, edim, etype, Physics>: CoordinateSystem, Physics {
	using CoordinateSystem::CoordinateSystem;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			double cos[3] { element.ecf.angle[gp][0], element.ecf.angle[gp][1], element.ecf.angle[gp][2] };
			double sin[3] { element.ecf.angle[gp][3], element.ecf.angle[gp][4], element.ecf.angle[gp][5] };
			double t[3][3] {
				{ cos[1] * cos[2]                           , cos[1] * sin[2]                           ,         -sin[1] },
				{ cos[2] * sin[0] * sin[1] - cos[0] * sin[2], cos[0] * cos[2] + sin[0] * sin[1] * sin[2], cos[1] * sin[0] },
				{ sin[0] * sin[2] + cos[0] * cos[2] * sin[1], cos[0] * sin[1] * sin[2] - cos[2] * sin[0], cos[0] * cos[1] }
			};

			double origin[9] = {
					element.ecf.conductivity[gp][0], element.ecf.conductivity[gp][1], element.ecf.conductivity[gp][2],
					element.ecf.conductivity[gp][3], element.ecf.conductivity[gp][4], element.ecf.conductivity[gp][5],
					element.ecf.conductivity[gp][6], element.ecf.conductivity[gp][7], element.ecf.conductivity[gp][8] };

			for (size_t i = 0; i < 3; ++i) {
				for (size_t j = 0; j < 3; ++j) {
					double a = t[0][i] * origin[0] + t[1][i] * origin[3] + t[2][i] * origin[6];
					double b = t[0][i] * origin[1] + t[1][i] * origin[4] + t[2][i] * origin[7];
					double c = t[0][i] * origin[2] + t[1][i] * origin[5] + t[2][i] * origin[8];
					element.conductivity[gp][3 * i + j] = a * t[0][j] + b * t[1][j] + c * t[2][j];
				}
			}
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cos0 = element.ecf.angle[gp][0];
			SIMD cos1 = element.ecf.angle[gp][1];
			SIMD cos2 = element.ecf.angle[gp][2];
			SIMD sin0 = element.ecf.angle[gp][3];
			SIMD sin1 = element.ecf.angle[gp][4];
			SIMD sin2 = element.ecf.angle[gp][5];

			SIMD t00 = cos1 * cos2;
			SIMD t01 = cos1 * sin2;
			SIMD t02 = -sin1;
			SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
			SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
			SIMD t12 = cos1 * sin0;
			SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
			SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
			SIMD t22 = cos0 * cos1;

			SIMD origin0 = element.ecf.conductivity[gp][9];
			SIMD origin1 = element.ecf.conductivity[gp][1];
			SIMD origin2 = element.ecf.conductivity[gp][2];
			SIMD origin3 = element.ecf.conductivity[gp][3];
			SIMD origin4 = element.ecf.conductivity[gp][4];
			SIMD origin5 = element.ecf.conductivity[gp][5];
			SIMD origin6 = element.ecf.conductivity[gp][6];
			SIMD origin7 = element.ecf.conductivity[gp][7];
			SIMD origin8 = element.ecf.conductivity[gp][8];

			SIMD a = t00 * origin0 + t10 * origin3 + t20 * origin6;
			SIMD b = t00 * origin1 + t10 * origin4 + t20 * origin7;
			SIMD c = t00 * origin2 + t10 * origin5 + t20 * origin8;
			element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

			a = t01 * origin0 + t11 * origin3 + t21 * origin6;
			b = t01 * origin1 + t11 * origin4 + t21 * origin7;
			c = t01 * origin2 + t11 * origin5 + t21 * origin8;
			element.conductivity[gp][3] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][4] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;

			a = t02 * origin0 + t12 * origin3 + t22 * origin6;
			b = t02 * origin1 + t12 * origin4 + t22 * origin7;
			c = t02 * origin2 + t12 * origin5 + t22 * origin8;
			element.conductivity[gp][6] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][7] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][8] = a * t02 + b * t12 + c * t22;
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_COORDINATESYSTEM_H_ */
