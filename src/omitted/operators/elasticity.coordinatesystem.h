
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_COORDINATESYSTEM_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_COORDINATESYSTEM_H_

#include <analysis/assembler/subkernel/operator.h>
#include "analysis/assembler/module/structuralmechanics.element.h"
#include "math/simd/simd.h"

#include <cmath>

namespace espreso {

struct ElasticityCoordinateSystem: ActionOperator {
	const char* name() const { return "ElasticityCoordinateSystem"; }

	ElasticityCoordinateSystem(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct ElasticityCoordinateSystemCartesian;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct ElasticityCoordinateSystemCylindric;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct ElasticityCoordinateSystemSpherical;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct ElasticityCoordinateSystemCopy;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct ElasticityCoordinateSystemApply;

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemCartesian<nodes, gps, 2, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	constexpr static double straightAngleRec = 1.0 / 180;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angle = element.ecf.center[gp][0];
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.cossin[gp][0][s] = std::cos(M_PI * angle[s] * straightAngleRec);
				element.cossin[gp][1][s] = std::sin(M_PI * angle[s] * straightAngleRec);
				element.cossin[gp][2][s] = std::cos(2 * M_PI * angle[s] * straightAngleRec);
				element.cossin[gp][3][s] = std::sin(2 * M_PI * angle[s] * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemCartesian<nodes, gps, 3, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	constexpr static double straightAngleRec = 1.0 / 180;

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angleX = element.ecf.center[gp][0];
			SIMD angleY = element.ecf.center[gp][1];
			SIMD angleZ = element.ecf.center[gp][2];
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.cossin[gp][ 0][s] = std::cos(M_PI * angleX[s] * straightAngleRec);
				element.cossin[gp][ 1][s] = std::cos(M_PI * angleY[s] * straightAngleRec);
				element.cossin[gp][ 2][s] = std::cos(M_PI * angleZ[s] * straightAngleRec);
				element.cossin[gp][ 3][s] = std::sin(M_PI * angleX[s] * straightAngleRec);
				element.cossin[gp][ 4][s] = std::sin(M_PI * angleY[s] * straightAngleRec);
				element.cossin[gp][ 5][s] = std::sin(M_PI * angleZ[s] * straightAngleRec);
				element.cossin[gp][ 6][s] = std::cos(2 * M_PI * angleX[s] * straightAngleRec);
				element.cossin[gp][ 7][s] = std::cos(2 * M_PI * angleY[s] * straightAngleRec);
				element.cossin[gp][ 8][s] = std::cos(2 * M_PI * angleZ[s] * straightAngleRec);
				element.cossin[gp][ 9][s] = std::sin(2 * M_PI * angleX[s] * straightAngleRec);
				element.cossin[gp][10][s] = std::sin(2 * M_PI * angleY[s] * straightAngleRec);
				element.cossin[gp][11][s] = std::sin(2 * M_PI * angleZ[s] * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemCylindric<nodes, gps, 2, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {

	ElasticityCoordinateSystemCylindric(size_t interval): ElasticityCoordinateSystem(interval)
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
				element.cossin[gp][2][s] = std::cos(2 * rot);
				element.cossin[gp][3][s] = std::sin(2 * rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemCylindric<nodes, gps, 3, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {

	ElasticityCoordinateSystemCylindric(size_t interval): ElasticityCoordinateSystem(interval)
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
				element.cossin[gp][ 0][s] = 1;
				element.cossin[gp][ 1][s] = 1;
				element.cossin[gp][ 2][s] = std::cos(rot);
				element.cossin[gp][ 3][s] = 0;
				element.cossin[gp][ 4][s] = 0;
				element.cossin[gp][ 5][s] = std::sin(rot);
				element.cossin[gp][ 6][s] = 1;
				element.cossin[gp][ 7][s] = 1;
				element.cossin[gp][ 8][s] = std::cos(2 * rot);
				element.cossin[gp][ 9][s] = 0;
				element.cossin[gp][10][s] = 0;
				element.cossin[gp][11][s] = std::sin(2 * rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemSpherical<nodes, gps, 3, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {

	ElasticityCoordinateSystemSpherical(size_t interval): ElasticityCoordinateSystem(interval)
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
				element.cossin[gp][ 0][s] = 1;
				element.cossin[gp][ 1][s] = std::cos(elevation);
				element.cossin[gp][ 2][s] = std::cos(azimut);
				element.cossin[gp][ 3][s] = 0;
				element.cossin[gp][ 4][s] = std::sin(elevation);
				element.cossin[gp][ 5][s] = std::sin(azimut);
				element.cossin[gp][ 6][s] = 1;
				element.cossin[gp][ 7][s] = std::cos(2 * elevation);
				element.cossin[gp][ 8][s] = std::cos(2 * azimut);
				element.cossin[gp][ 9][s] = 0;
				element.cossin[gp][10][s] = std::sin(2 * elevation);
				element.cossin[gp][11][s] = std::sin(2 * azimut);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ElasticityCoordinateSystemCopy<nodes, gps, 2, edim, StructuralMechanicsElementType::SYMMETRIC_PLANE, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		memcpy(element.elasticity, element.ecf.elasticity, sizeof(double) * SIMD::size * gps * 6);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ElasticityCoordinateSystemCopy<nodes, gps, 2, edim, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		memcpy(element.elasticity, element.ecf.elasticity, sizeof(double) * SIMD::size * gps * 10);
	}
};


template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemCopy<nodes, gps, 3, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		memcpy(element.elasticity, element.ecf.elasticity, sizeof(double) * SIMD::size * gps * 21);
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemApply<nodes, gps, 2, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	void simd(typename Physics::Element &element)
	{

	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemApply<nodes, gps, 3, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		SIMD C05 = load1(0.5);

		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD cosx  = element.cossin[gp][ 0];
			SIMD cosy  = element.cossin[gp][ 1];
			SIMD cosz  = element.cossin[gp][ 2];
			SIMD sinx  = element.cossin[gp][ 3];
			SIMD siny  = element.cossin[gp][ 4];
			SIMD sinz  = element.cossin[gp][ 5];
			SIMD cos2x = element.cossin[gp][ 6];
			SIMD cos2y = element.cossin[gp][ 7];
			SIMD cos2z = element.cossin[gp][ 8];
			SIMD sin2x = element.cossin[gp][ 9];
			SIMD sin2y = element.cossin[gp][10];
			SIMD sin2z = element.cossin[gp][11];

			SIMD T00 = cosx * cosx * cosz * cosz - sin2x * sin2z * cosy * C05 + cosy * cosy * sinx * sinx * sinz * sinz;
			SIMD T01 = cosz * cosz * sinx * sinx + sin2x * sin2z * cosy * C05 + cosx * cosx * cosy * cosy * sinz * sinz;
			SIMD T02 = siny * siny * sinz * sinz;
			SIMD T03 = sin2z * sinx * siny + sin2y * cosx * sinz * sinz;
			SIMD T04 = sin2z * cosx * siny - sin2y * sinx * sinz * sinz;
			SIMD T05 = sin2x * cosz * cosz + cos2x * sin2z * cosy - sin2x * cosy * cosy * sinz * sinz;
			SIMD T10 = cosx * cosx * sinz * sinz + sin2x * sin2z * cosy * C05 + cosy * cosy * cosz * cosz * sinx * sinx;
			SIMD T11 = sinx * sinx * sinz * sinz - sin2x * sin2z * cosy * C05 + cosx * cosx * cosy * cosy * cosz * cosz;
			SIMD T12 = cosz * cosz * siny * siny;
			SIMD T13 = -sin2z * sinx * siny + sin2y * cosx * cosz * cosz;
			SIMD T14 = -sin2z * cosx * siny - sin2y * cosz * cosz * sinx;
			SIMD T15 = sin2x * sinz * sinz - cos2x * sin2z * cosy - sin2x * cosy * cosy * cosz * cosz;
			SIMD T20 = sinx * sinx * siny * siny;
			SIMD T21 = cosx * cosx * siny * siny;
			SIMD T22 = cosy * cosy;
			SIMD T23 = -sin2y * cosx;
			SIMD T24 = sin2y * sinx;
			SIMD T25 = -sin2x * siny * siny;
			SIMD T30 = -C05 * sin2x * siny * sinz - sin2y * cosz * sinx * sinx * C05;
			SIMD T31 = sin2x * siny * sinz * C05 - sin2y * cosx * cosx * cosz * C05;
			SIMD T32 = sin2y * cosz * C05;
			SIMD T33 = cos2y * cosx * cosz - cosy * sinx * sinz;
			SIMD T34 = -cos2y * cosz * sinx - cosx * cosy * sinz;
			SIMD T35 = cos2x * siny * sinz + sin2x * sin2y * cosz * C05;
			SIMD T40 = sin2x * cosz * siny * C05 - sin2y * sinx * sinx * sinz * C05;
			SIMD T41 = -C05 * sin2x * cosz * siny - sin2y * cosx * cosx * sinz * C05;
			SIMD T42 = sin2y * sinz * C05;
			SIMD T43 = cos2y * cosx * sinz + cosy * cosz * sinx;
			SIMD T44 = -cos2y * sinx * sinz + cosx * cosy * cosz;
			SIMD T45 = -cos2x * cosz * siny + sin2x * sin2y * sinz * C05;
			SIMD T50 = -C05 * sin2z * cosx * cosx - cos2z * sin2x * cosy * C05 + sin2z * cosy * cosy * sinx * sinx * C05;
			SIMD T51 = -C05 * sin2z * sinx * sinx + cos2z * sin2x * cosy * C05 + sin2z * cosx * cosx * cosy * cosy * C05;
			SIMD T52 = sin2z * siny * siny * C05;
			SIMD T53 = cos2z * sinx * siny + sin2y * sin2z * cosx * C05;
			SIMD T54 = cos2z * cosx * siny - sin2y * sin2z * sinx * C05;
			SIMD T55 = -C05 * sin2x * sin2z + cos2x * cos2z * cosy - sin2x * sin2z * cosy * cosy * C05;

			SIMD e00 = element.ecf.elasticity[gp][ 0];
			SIMD e01 = element.ecf.elasticity[gp][ 1], e11 = element.ecf.elasticity[gp][ 6];
			SIMD e02 = element.ecf.elasticity[gp][ 2], e12 = element.ecf.elasticity[gp][ 7], e22 = element.ecf.elasticity[gp][11];
			SIMD e03 = element.ecf.elasticity[gp][ 3], e13 = element.ecf.elasticity[gp][ 8], e23 = element.ecf.elasticity[gp][12], e33 = element.ecf.elasticity[gp][15];
			SIMD e04 = element.ecf.elasticity[gp][ 4], e14 = element.ecf.elasticity[gp][ 9], e24 = element.ecf.elasticity[gp][13], e34 = element.ecf.elasticity[gp][16], e44 = element.ecf.elasticity[gp][18];
			SIMD e05 = element.ecf.elasticity[gp][ 5], e15 = element.ecf.elasticity[gp][10], e25 = element.ecf.elasticity[gp][14], e35 = element.ecf.elasticity[gp][17], e45 = element.ecf.elasticity[gp][19], e55 = element.ecf.elasticity[gp][20];

			SIMD a = T00 * e00 + T01 * e01 + T02 * e02 + T03 * e03 + T04 * e04 + T05 * e05;
			SIMD b = T00 * e01 + T01 * e11 + T02 * e12 + T03 * e13 + T04 * e14 + T05 * e15;
			SIMD c = T00 * e02 + T01 * e12 + T02 * e22 + T03 * e23 + T04 * e24 + T05 * e25;
			SIMD d = T00 * e03 + T01 * e13 + T02 * e23 + T03 * e33 + T04 * e34 + T05 * e35;
			SIMD e = T00 * e04 + T01 * e14 + T02 * e24 + T03 * e34 + T04 * e44 + T05 * e45;
			SIMD f = T00 * e05 + T01 * e15 + T02 * e25 + T03 * e35 + T04 * e45 + T05 * e55;
			element.elasticity[gp][ 0] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
			element.elasticity[gp][ 1] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][ 2] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][ 3] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][ 4] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][ 5] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T10 * e00 + T11 * e01 + T12 * e02 + T13 * e03 + T14 * e04 + T15 * e05;
			b = T10 * e01 + T11 * e11 + T12 * e12 + T13 * e13 + T14 * e14 + T15 * e15;
			c = T10 * e02 + T11 * e12 + T12 * e22 + T13 * e23 + T14 * e24 + T15 * e25;
			d = T10 * e03 + T11 * e13 + T12 * e23 + T13 * e33 + T14 * e34 + T15 * e35;
			e = T10 * e04 + T11 * e14 + T12 * e24 + T13 * e34 + T14 * e44 + T15 * e45;
			f = T10 * e05 + T11 * e15 + T12 * e25 + T13 * e35 + T14 * e45 + T15 * e55;
			element.elasticity[gp][ 6] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][ 7] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][ 8] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][ 9] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][10] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T20 * e00 + T21 * e01 + T22 * e02 + T23 * e03 + T24 * e04 + T25 * e05;
			b = T20 * e01 + T21 * e11 + T22 * e12 + T23 * e13 + T24 * e14 + T25 * e15;
			c = T20 * e02 + T21 * e12 + T22 * e22 + T23 * e23 + T24 * e24 + T25 * e25;
			d = T20 * e03 + T21 * e13 + T22 * e23 + T23 * e33 + T24 * e34 + T25 * e35;
			e = T20 * e04 + T21 * e14 + T22 * e24 + T23 * e34 + T24 * e44 + T25 * e45;
			f = T20 * e05 + T21 * e15 + T22 * e25 + T23 * e35 + T24 * e45 + T25 * e55;
			element.elasticity[gp][11] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][12] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][13] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][14] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T30 * e00 + T31 * e01 + T32 * e02 + T33 * e03 + T34 * e04 + T35 * e05;
			b = T30 * e01 + T31 * e11 + T32 * e12 + T33 * e13 + T34 * e14 + T35 * e15;
			c = T30 * e02 + T31 * e12 + T32 * e22 + T33 * e23 + T34 * e24 + T35 * e25;
			d = T30 * e03 + T31 * e13 + T32 * e23 + T33 * e33 + T34 * e34 + T35 * e35;
			e = T30 * e04 + T31 * e14 + T32 * e24 + T33 * e34 + T34 * e44 + T35 * e45;
			f = T30 * e05 + T31 * e15 + T32 * e25 + T33 * e35 + T34 * e45 + T35 * e55;
			element.elasticity[gp][15] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][16] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][17] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T40 * e00 + T41 * e01 + T42 * e02 + T43 * e03 + T44 * e04 + T45 * e05;
			b = T40 * e01 + T41 * e11 + T42 * e12 + T43 * e13 + T44 * e14 + T45 * e15;
			c = T40 * e02 + T41 * e12 + T42 * e22 + T43 * e23 + T44 * e24 + T45 * e25;
			d = T40 * e03 + T41 * e13 + T42 * e23 + T43 * e33 + T44 * e34 + T45 * e35;
			e = T40 * e04 + T41 * e14 + T42 * e24 + T43 * e34 + T44 * e44 + T45 * e45;
			f = T40 * e05 + T41 * e15 + T42 * e25 + T43 * e35 + T44 * e45 + T45 * e55;
			element.elasticity[gp][18] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][19] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T50 * e00 + T51 * e01 + T52 * e02 + T53 * e03 + T54 * e04 + T55 * e05;
			b = T50 * e01 + T51 * e11 + T52 * e12 + T53 * e13 + T54 * e14 + T55 * e15;
			c = T50 * e02 + T51 * e12 + T52 * e22 + T53 * e23 + T54 * e24 + T55 * e25;
			d = T50 * e03 + T51 * e13 + T52 * e23 + T53 * e33 + T54 * e34 + T55 * e35;
			e = T50 * e04 + T51 * e14 + T52 * e24 + T53 * e34 + T54 * e44 + T55 * e45;
			f = T50 * e05 + T51 * e15 + T52 * e25 + T53 * e35 + T54 * e45 + T55 * e55;
			element.elasticity[gp][20] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_COORDINATESYSTEM_H_ */
