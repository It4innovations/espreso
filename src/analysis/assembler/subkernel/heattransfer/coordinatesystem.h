
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_COORDINATESYSTEM_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_COORDINATESYSTEM_H_

#include "subkernels.h"
#include "conductivity.h"
#include <cmath>

namespace espreso {

struct HeatTransferCoordinateSystem: SubKernel {
	const char* name() const { return "HeatTransferCoordinateSystemKernel"; }

	const CoordinateSystemConfiguration *configuration;
	CoordinateSystemConfiguration::TYPE type;
	bool rotated;

	HeatTransferCoordinateSystem()
	: configuration(nullptr), type(CoordinateSystemConfiguration::TYPE::CARTESIAN), rotated(false)
	{
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE;
	}

	void activate(const CoordinateSystemConfiguration &configuration, bool isconst, bool rotated)
	{
		this->configuration = &configuration;
		this->type = configuration.type;
		this->isconst = isconst;
		this->rotated = rotated;
		this->isactive = 1;
		switch (this->type) {
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: this->isconst = 0; break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:   this->isconst = 0; break;
		}
	}
};

template <size_t gps, size_t ndim, class Physics> struct HeatTransferCoordinateSystemCartesian;
template <size_t gps, size_t ndim, class Physics> struct HeatTransferCoordinateSystemCylindric;

template <size_t gps, class Physics>
struct HeatTransferCoordinateSystemCartesian<gps, 2, Physics> {

	constexpr static double straightAngleRec = 1.0 / 180;

	static void simd(typename Physics::Element &element)
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

template <size_t gps, class Physics>
struct HeatTransferCoordinateSystemCartesian<gps, 3, Physics> {

	constexpr static double straightAngleRec = 1.0 / 180;

	static void simd(typename Physics::Element &element)
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

template <size_t gps, class Physics>
struct HeatTransferCoordinateSystemCylindric<gps, 2, Physics> {

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
				element.cossin[gp][0][s] = std::cos(rot);
				element.cossin[gp][1][s] = std::sin(rot);
			}
		}
	}
};

template <size_t gps, class Physics>
struct HeatTransferCoordinateSystemCylindric<gps, 3, Physics> {

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

template <size_t gps, class Physics>
struct HeatTransferCoordinateSystemSpherical {

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

template <size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model, class Physics> struct HeatTransferCoordinateSystemKernel: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element) {}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 2, ThermalConductivityConfiguration::MODEL::ISOTROPIC, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
			element.conductivity[gps][1] = element.ecf.conductivity[gps][0];
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 3, ThermalConductivityConfiguration::MODEL::ISOTROPIC, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
			element.conductivity[gps][1] = element.ecf.conductivity[gps][0];
			element.conductivity[gps][2] = element.ecf.conductivity[gps][0];
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 2, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
			element.conductivity[gps][1] = element.ecf.conductivity[gps][1];
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 3, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::DIAGONAL, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
			element.conductivity[gps][1] = element.ecf.conductivity[gps][1];
			element.conductivity[gps][2] = element.ecf.conductivity[gps][2];
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 2, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::SYMMETRIC, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		switch (type) { // always rotated coordinate system
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:   HeatTransferCoordinateSystemCartesian<gps, 2, Physics>::simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: HeatTransferCoordinateSystemCylindric<gps, 2, Physics>::simd(element); break;
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			SIMD c00 = element.ecf.conductivity[gp][0];
			SIMD c11 = element.ecf.conductivity[gp][1];
			SIMD cos = element.cossin[gp][0];
			SIMD sin = element.cossin[gp][1];

			element.conductivity[gp][0] = (cos * c00) * cos + (sin * c11) * sin;
			element.conductivity[gp][1] = (cos * c00) * sin - (sin * c11) * cos;
			element.conductivity[gp][2] = (sin * c00) * sin + (cos * c11) * cos;
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 3, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::SYMMETRIC, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		switch (type) { // always rotated coordinate system
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:   HeatTransferCoordinateSystemCartesian<gps, 3, Physics>::simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: HeatTransferCoordinateSystemCylindric<gps, 3, Physics>::simd(element); break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:   HeatTransferCoordinateSystemSpherical<gps, Physics>::simd(element); break;
		}
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
			SIMD c11 = element.ecf.conductivity[gp][1];
			SIMD c22 = element.ecf.conductivity[gp][2];

			SIMD a = t00 * c00, b = t10 * c11, c = t20 * c22;
			element.conductivity[gp][0] = a * t00 + b * t10 + c * t20;
			element.conductivity[gp][1] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][2] = a * t02 + b * t12 + c * t22;

			a = t01 * c00, b = t11 * c11, c = t21 * c22;
			element.conductivity[gp][3] = a * t01 + b * t11 + c * t21;
			element.conductivity[gp][4] = a * t02 + b * t12 + c * t22;

			a = t02 * c00, b = t12 * c11, c = t22 * c22;
			element.conductivity[gp][5] = a * t02 + b * t12 + c * t22;
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 2, ThermalConductivityConfiguration::MODEL::SYMMETRIC, ThermalConductivityConfiguration::MODEL::SYMMETRIC, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		if (rotated) {
			switch (type) { // always rotated coordinate system
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:   HeatTransferCoordinateSystemCartesian<gps, 2, Physics>::simd(element); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: HeatTransferCoordinateSystemCylindric<gps, 2, Physics>::simd(element); break;
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD c00 = element.ecf.conductivity[gp][0];
				SIMD c01 = element.ecf.conductivity[gp][1], c11 = element.ecf.conductivity[gp][2];
				SIMD cos = element.cossin[gp][0];
				SIMD sin = element.cossin[gp][1];

				element.conductivity[gp][0] = (cos * c00 - sin * c01) * cos - (cos * c01 - sin * c11) * sin;
				element.conductivity[gp][1] = (cos * c00 - sin * c01) * sin + (cos * c01 - sin * c11) * cos;
				element.conductivity[gp][2] = (sin * c00 + cos * c01) * sin + (sin * c01 + cos * c11) * cos;
			}
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
				element.conductivity[gps][1] = element.ecf.conductivity[gps][1];
				element.conductivity[gps][2] = element.ecf.conductivity[gps][2];
			}
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 3, ThermalConductivityConfiguration::MODEL::SYMMETRIC, ThermalConductivityConfiguration::MODEL::SYMMETRIC, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		if (rotated) {
			switch (type) { // always rotated coordinate system
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:   HeatTransferCoordinateSystemCartesian<gps, 3, Physics>::simd(element); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: HeatTransferCoordinateSystemCylindric<gps, 3, Physics>::simd(element); break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:   HeatTransferCoordinateSystemSpherical<gps, Physics>::simd(element); break;
			}
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
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
				element.conductivity[gps][1] = element.ecf.conductivity[gps][1];
				element.conductivity[gps][2] = element.ecf.conductivity[gps][2];
				element.conductivity[gps][3] = element.ecf.conductivity[gps][3];
				element.conductivity[gps][4] = element.ecf.conductivity[gps][4];
				element.conductivity[gps][5] = element.ecf.conductivity[gps][5];
			}
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 2, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		if (rotated) {
			switch (type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:   HeatTransferCoordinateSystemCartesian<gps, 2, Physics>::simd(element); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: HeatTransferCoordinateSystemCylindric<gps, 2, Physics>::simd(element); break;
			}

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
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
				element.conductivity[gps][1] = element.ecf.conductivity[gps][1];
				element.conductivity[gps][2] = element.ecf.conductivity[gps][2];
				element.conductivity[gps][3] = element.ecf.conductivity[gps][3];
			}
		}
	}
};

template <size_t gps, class Physics> struct HeatTransferCoordinateSystemKernel<gps, 3, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, Physics>: HeatTransferCoordinateSystem, Physics {
	HeatTransferCoordinateSystemKernel(const HeatTransferCoordinateSystem &base): HeatTransferCoordinateSystem(base) {}

	void simd(typename Physics::Element &element)
	{
		if (rotated) {
			switch (type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:   HeatTransferCoordinateSystemCartesian<gps, 3, Physics>::simd(element); break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: HeatTransferCoordinateSystemCylindric<gps, 3, Physics>::simd(element); break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:   HeatTransferCoordinateSystemSpherical<gps, Physics>::simd(element); break;
			}

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
		} else {
			for (size_t gp = 0; gp < gps; ++gp) {
				element.conductivity[gps][0] = element.ecf.conductivity[gps][0];
				element.conductivity[gps][1] = element.ecf.conductivity[gps][1];
				element.conductivity[gps][2] = element.ecf.conductivity[gps][2];
				element.conductivity[gps][3] = element.ecf.conductivity[gps][3];
				element.conductivity[gps][4] = element.ecf.conductivity[gps][4];
				element.conductivity[gps][5] = element.ecf.conductivity[gps][5];
				element.conductivity[gps][6] = element.ecf.conductivity[gps][6];
				element.conductivity[gps][7] = element.ecf.conductivity[gps][7];
				element.conductivity[gps][8] = element.ecf.conductivity[gps][8];
				element.conductivity[gps][9] = element.ecf.conductivity[gps][9];
			}
		}
	}
};
}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_COORDINATESYSTEM_H_ */
