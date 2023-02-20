
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_COORDINATESYSTEM_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_COORDINATESYSTEM_H_

#include "analysis/assembler/operator.h"
#include "math/simd/simd.h"

#include <cmath>

namespace espreso {

struct ElasticityCoordinateSystem: ActionOperator {
	ElasticityCoordinateSystem(size_t interval)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemCopy<nodes, gps, 2, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		memcpy(element.elasticity, element.ecf.elasticity, sizeof(double) * SIMD::size * gps * 16);
	}
};

template <size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct ElasticityCoordinateSystemCopy<nodes, gps, 3, edim, etype, Physics>: ElasticityCoordinateSystem, Physics {
	using ElasticityCoordinateSystem::ElasticityCoordinateSystem;

	void simd(typename Physics::Element &element)
	{
		memcpy(element.elasticity, element.ecf.elasticity, sizeof(double) * SIMD::size * gps * 36);
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
			SIMD e01 = element.ecf.elasticity[gp][ 1];
			SIMD e02 = element.ecf.elasticity[gp][ 2];
			SIMD e03 = element.ecf.elasticity[gp][ 3];
			SIMD e04 = element.ecf.elasticity[gp][ 4];
			SIMD e05 = element.ecf.elasticity[gp][ 5];
			SIMD e10 = element.ecf.elasticity[gp][ 6];
			SIMD e11 = element.ecf.elasticity[gp][ 7];
			SIMD e12 = element.ecf.elasticity[gp][ 8];
			SIMD e13 = element.ecf.elasticity[gp][ 9];
			SIMD e14 = element.ecf.elasticity[gp][10];
			SIMD e15 = element.ecf.elasticity[gp][11];
			SIMD e20 = element.ecf.elasticity[gp][12];
			SIMD e21 = element.ecf.elasticity[gp][13];
			SIMD e22 = element.ecf.elasticity[gp][14];
			SIMD e23 = element.ecf.elasticity[gp][15];
			SIMD e24 = element.ecf.elasticity[gp][16];
			SIMD e25 = element.ecf.elasticity[gp][17];
			SIMD e30 = element.ecf.elasticity[gp][18];
			SIMD e31 = element.ecf.elasticity[gp][19];
			SIMD e32 = element.ecf.elasticity[gp][20];
			SIMD e33 = element.ecf.elasticity[gp][21];
			SIMD e34 = element.ecf.elasticity[gp][22];
			SIMD e35 = element.ecf.elasticity[gp][23];
			SIMD e40 = element.ecf.elasticity[gp][24];
			SIMD e41 = element.ecf.elasticity[gp][25];
			SIMD e42 = element.ecf.elasticity[gp][26];
			SIMD e43 = element.ecf.elasticity[gp][27];
			SIMD e44 = element.ecf.elasticity[gp][28];
			SIMD e45 = element.ecf.elasticity[gp][29];
			SIMD e50 = element.ecf.elasticity[gp][30];
			SIMD e51 = element.ecf.elasticity[gp][31];
			SIMD e52 = element.ecf.elasticity[gp][32];
			SIMD e53 = element.ecf.elasticity[gp][33];
			SIMD e54 = element.ecf.elasticity[gp][34];
			SIMD e55 = element.ecf.elasticity[gp][35];

			SIMD a = T00 * e00 + T01 * e10 + T02 * e20 + T03 * e30 + T04 * e40 + T05 * e50;
			SIMD b = T00 * e01 + T01 * e11 + T02 * e21 + T03 * e31 + T04 * e41 + T05 * e51;
			SIMD c = T00 * e02 + T01 * e12 + T02 * e22 + T03 * e32 + T04 * e42 + T05 * e52;
			SIMD d = T00 * e03 + T01 * e13 + T02 * e23 + T03 * e33 + T04 * e43 + T05 * e53;
			SIMD e = T00 * e04 + T01 * e14 + T02 * e24 + T03 * e34 + T04 * e44 + T05 * e54;
			SIMD f = T00 * e05 + T01 * e15 + T02 * e25 + T03 * e35 + T04 * e45 + T05 * e55;
			element.elasticity[gp][ 0] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
			element.elasticity[gp][ 1] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][ 2] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][ 3] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][ 4] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][ 5] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T10 * e00 + T11 * e10 + T12 * e20 + T13 * e30 + T14 * e40 + T15 * e50;
			b = T10 * e01 + T11 * e11 + T12 * e21 + T13 * e31 + T14 * e41 + T15 * e51;
			c = T10 * e02 + T11 * e12 + T12 * e22 + T13 * e32 + T14 * e42 + T15 * e52;
			d = T10 * e03 + T11 * e13 + T12 * e23 + T13 * e33 + T14 * e43 + T15 * e53;
			e = T10 * e04 + T11 * e14 + T12 * e24 + T13 * e34 + T14 * e44 + T15 * e54;
			f = T10 * e05 + T11 * e15 + T12 * e25 + T13 * e35 + T14 * e45 + T15 * e55;
			element.elasticity[gp][ 6] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
			element.elasticity[gp][ 7] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][ 8] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][ 9] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][10] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][11] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T20 * e00 + T21 * e10 + T22 * e20 + T23 * e30 + T24 * e40 + T25 * e50;
			b = T20 * e01 + T21 * e11 + T22 * e21 + T23 * e31 + T24 * e41 + T25 * e51;
			c = T20 * e02 + T21 * e12 + T22 * e22 + T23 * e32 + T24 * e42 + T25 * e52;
			d = T20 * e03 + T21 * e13 + T22 * e23 + T23 * e33 + T24 * e43 + T25 * e53;
			e = T20 * e04 + T21 * e14 + T22 * e24 + T23 * e34 + T24 * e44 + T25 * e54;
			f = T20 * e05 + T21 * e15 + T22 * e25 + T23 * e35 + T24 * e45 + T25 * e55;
			element.elasticity[gp][12] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
			element.elasticity[gp][13] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][14] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][15] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][16] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][17] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T30 * e00 + T31 * e10 + T32 * e20 + T33 * e30 + T34 * e40 + T35 * e50;
			b = T30 * e01 + T31 * e11 + T32 * e21 + T33 * e31 + T34 * e41 + T35 * e51;
			c = T30 * e02 + T31 * e12 + T32 * e22 + T33 * e32 + T34 * e42 + T35 * e52;
			d = T30 * e03 + T31 * e13 + T32 * e23 + T33 * e33 + T34 * e43 + T35 * e53;
			e = T30 * e04 + T31 * e14 + T32 * e24 + T33 * e34 + T34 * e44 + T35 * e54;
			f = T30 * e05 + T31 * e15 + T32 * e25 + T33 * e35 + T34 * e45 + T35 * e55;
			element.elasticity[gp][18] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
			element.elasticity[gp][19] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][20] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][21] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][22] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][23] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T40 * e00 + T41 * e10 + T42 * e20 + T43 * e30 + T44 * e40 + T45 * e50;
			b = T40 * e01 + T41 * e11 + T42 * e21 + T43 * e31 + T44 * e41 + T45 * e51;
			c = T40 * e02 + T41 * e12 + T42 * e22 + T43 * e32 + T44 * e42 + T45 * e52;
			d = T40 * e03 + T41 * e13 + T42 * e23 + T43 * e33 + T44 * e43 + T45 * e53;
			e = T40 * e04 + T41 * e14 + T42 * e24 + T43 * e34 + T44 * e44 + T45 * e54;
			f = T40 * e05 + T41 * e15 + T42 * e25 + T43 * e35 + T44 * e45 + T45 * e55;
			element.elasticity[gp][24] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
			element.elasticity[gp][25] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][26] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][27] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][28] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][29] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;

			a = T50 * e00 + T51 * e10 + T52 * e20 + T53 * e30 + T54 * e40 + T55 * e50;
			b = T50 * e01 + T51 * e11 + T52 * e21 + T53 * e31 + T54 * e41 + T55 * e51;
			c = T50 * e02 + T51 * e12 + T52 * e22 + T53 * e32 + T54 * e42 + T55 * e52;
			d = T50 * e03 + T51 * e13 + T52 * e23 + T53 * e33 + T54 * e43 + T55 * e53;
			e = T50 * e04 + T51 * e14 + T52 * e24 + T53 * e34 + T54 * e44 + T55 * e54;
			f = T50 * e05 + T51 * e15 + T52 * e25 + T53 * e35 + T54 * e45 + T55 * e55;
			element.elasticity[gp][30] = a * T00 + b * T01 + c * T02 + d * T03 + e * T04 + f * T05;
			element.elasticity[gp][31] = a * T10 + b * T11 + c * T12 + d * T13 + e * T14 + f * T15;
			element.elasticity[gp][32] = a * T20 + b * T21 + c * T22 + d * T23 + e * T24 + f * T25;
			element.elasticity[gp][33] = a * T30 + b * T31 + c * T32 + d * T33 + e * T34 + f * T35;
			element.elasticity[gp][34] = a * T40 + b * T41 + c * T42 + d * T43 + e * T44 + f * T45;
			element.elasticity[gp][35] = a * T50 + b * T51 + c * T52 + d * T53 + e * T54 + f * T55;
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ELASTICITY_COORDINATESYSTEM_H_ */
