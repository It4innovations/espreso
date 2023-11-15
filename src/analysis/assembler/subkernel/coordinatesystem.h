
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATESYSTEM_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATESYSTEM_H_

#include <cmath>

namespace espreso {

template <size_t gps, size_t ndim, size_t multiplicity, class Physics> struct CoordinateSystemCartesian;
template <size_t gps, size_t ndim, size_t multiplicity, class Physics> struct CoordinateSystemCylindric;

template <size_t gps, size_t multiplicity, class Physics>
struct CoordinateSystemCartesian<gps, 2, multiplicity, Physics> {

	constexpr static double straightAngleRec = 1.0 / 180;

	static void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angle = element.ecf.center[gp][0];
			for (size_t s = 0; s < SIMD::size; ++s) {
				for (size_t m = 1; m <= multiplicity; ++m) {
					element.cossin[gp][(m - 1) * 2 + 0][s] = std::cos(m * M_PI * angle[s] * straightAngleRec);
					element.cossin[gp][(m - 1) * 2 + 1][s] = std::sin(m * M_PI * angle[s] * straightAngleRec);
				}
			}
		}
	}
};

template <size_t gps, size_t multiplicity, class Physics>
struct CoordinateSystemCartesian<gps, 3, multiplicity, Physics> {

	constexpr static double straightAngleRec = 1.0 / 180;

	static void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD angleX = element.ecf.center[gp][0];
			SIMD angleY = element.ecf.center[gp][1];
			SIMD angleZ = element.ecf.center[gp][2];
			for (size_t s = 0; s < SIMD::size; ++s) {
				for (size_t m = 1; m <= multiplicity; ++m) {
					element.cossin[gp][(m - 1) * 6 + 0][s] = std::cos(m * M_PI * angleX[s] * straightAngleRec);
					element.cossin[gp][(m - 1) * 6 + 1][s] = std::cos(m * M_PI * angleY[s] * straightAngleRec);
					element.cossin[gp][(m - 1) * 6 + 2][s] = std::cos(m * M_PI * angleZ[s] * straightAngleRec);
					element.cossin[gp][(m - 1) * 6 + 3][s] = std::sin(m * M_PI * angleX[s] * straightAngleRec);
					element.cossin[gp][(m - 1) * 6 + 4][s] = std::sin(m * M_PI * angleY[s] * straightAngleRec);
					element.cossin[gp][(m - 1) * 6 + 5][s] = std::sin(m * M_PI * angleZ[s] * straightAngleRec);
				}
			}
		}
	}
};

template <size_t gps, size_t multiplicity, class Physics>
struct CoordinateSystemCylindric<gps, 2, multiplicity, Physics> {

	static void simd(typename Physics::Element &element)
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
				for (size_t m = 1; m <= multiplicity; ++m) {
					element.cossin[gp][(m - 1) * 2 + 0][s] = std::cos(m * rot);
					element.cossin[gp][(m - 1) * 2 + 1][s] = std::sin(m * rot);
				}
			}
		}
	}
};

template <size_t gps, size_t multiplicity, class Physics>
struct CoordinateSystemCylindric<gps, 3, multiplicity, Physics> {

	static void simd(typename Physics::Element &element)
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
				for (size_t m = 1; m <= multiplicity; ++m) {
					element.cossin[gp][(m - 1) * 6 + 0][s] = 1;
					element.cossin[gp][(m - 1) * 6 + 1][s] = 1;
					element.cossin[gp][(m - 1) * 6 + 2][s] = std::cos(m * rot);
					element.cossin[gp][(m - 1) * 6 + 3][s] = 0;
					element.cossin[gp][(m - 1) * 6 + 4][s] = 0;
					element.cossin[gp][(m - 1) * 6 + 5][s] = std::sin(m * rot);
				}
			}
		}
	}
};

template <size_t gps, size_t multiplicity, class Physics>
struct CoordinateSystemSpherical {

	static void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD x = element.gpcoords[gp][0] - element.ecf.center[gp][0];
			SIMD y = element.gpcoords[gp][1] - element.ecf.center[gp][1];
			SIMD z = element.gpcoords[gp][2] - element.ecf.center[gp][2];
			for (size_t s = 0; s < SIMD::size; ++s) {
				double azimut = std::atan2(y[s], x[s]);
				double r = std::sqrt(x[s] * x[s] + y[s] * y[s] + z[s] * z[s]);
				double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z[s] * z[s] + x[s] * x[s]), y[s]);
				for (size_t m = 1; m <= multiplicity; ++m) {
					element.cossin[gp][(m - 1) * 6 + 0][s] = 1;
					element.cossin[gp][(m - 1) * 6 + 1][s] = std::cos(m * elevation);
					element.cossin[gp][(m - 1) * 6 + 2][s] = std::cos(m * azimut);
					element.cossin[gp][(m - 1) * 6 + 3][s] = 0;
					element.cossin[gp][(m - 1) * 6 + 4][s] = std::sin(m * elevation);
					element.cossin[gp][(m - 1) * 6 + 5][s] = std::sin(m * azimut);
				}
			}
		}
	}
};

template <size_t gps, size_t ndim, size_t multiplicity, class Physics> struct CoordinateSystem;

template <size_t gps, size_t multiplicity, class Physics> struct CoordinateSystem<gps, 2, multiplicity, Physics> {

	static void simd(typename Physics::Element &element, CoordinateSystemConfiguration::TYPE type, int isconst)
	{
		switch (type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:   if (!isconst) CoordinateSystemCartesian<gps, 2, multiplicity, Physics>::simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:               CoordinateSystemCylindric<gps, 2, multiplicity, Physics>::simd(element); break;
		}
	}
};

template <size_t gps, size_t multiplicity, class Physics> struct CoordinateSystem<gps, 3, multiplicity, Physics> {

	static void simd(typename Physics::Element &element, CoordinateSystemConfiguration::TYPE type, int isconst)
	{
		switch (type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:   if (!isconst) CoordinateSystemCartesian<gps, 3, multiplicity, Physics>::simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:               CoordinateSystemCylindric<gps, 3, multiplicity, Physics>::simd(element); break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:                 CoordinateSystemSpherical<gps,    multiplicity, Physics>::simd(element); break;
		}
	}
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_COORDINATESYSTEM_H_ */
