
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

namespace espreso {

struct Conductivity: ActionOperator {
	Conductivity(int interval) {}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct ConductivityDiagonal;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct ConductivitySymmetric;

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityDiagonal<nodes, gps, 2, edim, Physics>: Conductivity, Physics {
	using Conductivity::Conductivity;

	// [ KXX, KYY ]     [ KXX,  0  ]
	// [ ---  --- ]  -> [  0   KYY ]
	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			element.ecf.conductivity[4 * gpindex + 3] = element.ecf.conductivity[4 * gpindex + 1];
			element.ecf.conductivity[4 * gpindex + 1] = element.ecf.conductivity[4 * gpindex + 2] = 0;
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			store(element.ecf.conductivity + (4 * gpindex + 3) * SIMD::size, load(element.ecf.conductivity + (4 * gpindex + 1) * SIMD::size));
			store(element.ecf.conductivity + (4 * gpindex + 1) * SIMD::size, zeros());
			store(element.ecf.conductivity + (4 * gpindex + 2) * SIMD::size, zeros());
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityDiagonal<nodes, gps, 3, edim, Physics>: Conductivity, Physics {
	using Conductivity::Conductivity;

	// [ KXX, KYY, KZZ ]     [ KXX,  0  , 0  ]
	// [ ---  ---  --- ]  -> [  0   KYY   0  ]
	// [ ---  ---  --- ]     [  0    0   KZZ ]
	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			element.ecf.conductivity[9 * gpindex + 4] = element.ecf.conductivity[9 * gpindex + 1];
			element.ecf.conductivity[9 * gpindex + 8] = element.ecf.conductivity[9 * gpindex + 2];
			element.ecf.conductivity[9 * gpindex + 1] = element.ecf.conductivity[9 * gpindex + 2] = 0;
			element.ecf.conductivity[9 * gpindex + 3] = element.ecf.conductivity[9 * gpindex + 5] = 0;
			element.ecf.conductivity[9 * gpindex + 6] = element.ecf.conductivity[9 * gpindex + 7] = 0;
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			store(element.ecf.conductivity + (9 * gpindex + 4) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 1) * SIMD::size));
			store(element.ecf.conductivity + (9 * gpindex + 8) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 2) * SIMD::size));
			store(element.ecf.conductivity + (9 * gpindex + 1) * SIMD::size, zeros());
			store(element.ecf.conductivity + (9 * gpindex + 2) * SIMD::size, zeros());
			store(element.ecf.conductivity + (9 * gpindex + 3) * SIMD::size, zeros());
			store(element.ecf.conductivity + (9 * gpindex + 5) * SIMD::size, zeros());
			store(element.ecf.conductivity + (9 * gpindex + 6) * SIMD::size, zeros());
			store(element.ecf.conductivity + (9 * gpindex + 7) * SIMD::size, zeros());
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivitySymmetric<nodes, gps, 2, edim, Physics>: Conductivity, Physics {
	using Conductivity::Conductivity;

	// [ KXX, KXY ]     [ KXX, KXY ]
	// [ KYY  --- ]  -> [ KXY  KYY ]
	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			element.ecf.conductivity[4 * gpindex + 3] = element.ecf.conductivity[4 * gpindex + 2];
			element.ecf.conductivity[4 * gpindex + 2] = element.ecf.conductivity[4 * gpindex + 1];
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			store(element.ecf.conductivity + (4 * gpindex + 3) * SIMD::size, load(element.ecf.conductivity + (4 * gpindex + 2) * SIMD::size));
			store(element.ecf.conductivity + (4 * gpindex + 2) * SIMD::size, load(element.ecf.conductivity + (4 * gpindex + 1) * SIMD::size));
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivitySymmetric<nodes, gps, 3, edim, Physics>: Conductivity, Physics {
	using Conductivity::Conductivity;

	// [ KXX, KXY, KXZ ]     [ KXX, KXY ,KXZ ]  [ 0, 1 ,2 ]
	// [ KYY  KYZ  KZZ ]  -> [ KXY  KYY  KYZ ]  [ 3  4  5 ]
	// [ ---  ---  --- ]     [ KXZ  KYZ  KZZ ]  [ 6  7  8 ]
	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			element.ecf.conductivity[9 * gpindex + 8] = element.ecf.conductivity[9 * gpindex + 5];
			element.ecf.conductivity[9 * gpindex + 7] = element.ecf.conductivity[9 * gpindex + 4];
			element.ecf.conductivity[9 * gpindex + 6] = element.ecf.conductivity[9 * gpindex + 2];

			element.ecf.conductivity[9 * gpindex + 5] = element.ecf.conductivity[9 * gpindex + 4];
			element.ecf.conductivity[9 * gpindex + 4] = element.ecf.conductivity[9 * gpindex + 3];
			element.ecf.conductivity[9 * gpindex + 3] = element.ecf.conductivity[9 * gpindex + 1];
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			store(element.ecf.conductivity + (9 * gpindex + 8) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 5) * SIMD::size));
			store(element.ecf.conductivity + (9 * gpindex + 7) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 4) * SIMD::size));
			store(element.ecf.conductivity + (9 * gpindex + 6) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 2) * SIMD::size));
			store(element.ecf.conductivity + (9 * gpindex + 5) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 4) * SIMD::size));
			store(element.ecf.conductivity + (9 * gpindex + 4) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 3) * SIMD::size));
			store(element.ecf.conductivity + (9 * gpindex + 3) * SIMD::size, load(element.ecf.conductivity + (9 * gpindex + 1) * SIMD::size));

		}
	}
};

struct ConductivityRotation: ActionOperator {
	ConductivityRotation(int interval) {}
};


template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct ConductivityRotationCartesian;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct ConductivityRotationCylindric;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct ConductivityRotationSpherical;
template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct ConductivityRotationApply;

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityRotationCartesian<nodes, gps, 2, edim, Physics>: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	constexpr static double straightAngleRec = 1.0 / 180;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double angle = element.ecf.angle[2 * gpindex + 0];
			element.ecf.angle[2 * gpindex + 0] = std::cos(M_PI * angle * straightAngleRec);
			element.ecf.angle[2 * gpindex + 1] = std::sin(M_PI * angle * straightAngleRec);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				double angle = element.ecf.angle[(2 * gpindex + 0) * SIMD::size + s];
				element.ecf.angle[(2 * gpindex + 0) * SIMD::size + s] = std::cos(M_PI * angle * straightAngleRec);
				element.ecf.angle[(2 * gpindex + 1) * SIMD::size + s] = std::sin(M_PI * angle * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityRotationCartesian<nodes, gps, 3, edim, Physics>: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	constexpr static double straightAngleRec = 1.0 / 180;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double angleX = element.ecf.angle[6 * gpindex + 0];
			double angleY = element.ecf.angle[6 * gpindex + 1];
			double angleZ = element.ecf.angle[6 * gpindex + 2];
			element.ecf.angle[6 * gpindex + 0] = std::cos(M_PI * angleX * straightAngleRec);
			element.ecf.angle[6 * gpindex + 1] = std::cos(M_PI * angleY * straightAngleRec);
			element.ecf.angle[6 * gpindex + 2] = std::cos(M_PI * angleZ * straightAngleRec);
			element.ecf.angle[6 * gpindex + 3] = std::sin(M_PI * angleX * straightAngleRec);
			element.ecf.angle[6 * gpindex + 4] = std::sin(M_PI * angleY * straightAngleRec);
			element.ecf.angle[6 * gpindex + 5] = std::sin(M_PI * angleZ * straightAngleRec);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				double angleX = element.ecf.angle[(6 * gpindex + 0) * SIMD::size + s];
				double angleY = element.ecf.angle[(6 * gpindex + 1) * SIMD::size + s];
				double angleZ = element.ecf.angle[(6 * gpindex + 2) * SIMD::size + s];
				element.ecf.angle[(6 * gpindex + 0) * SIMD::size + s] = std::cos(M_PI * angleX * straightAngleRec);
				element.ecf.angle[(6 * gpindex + 1) * SIMD::size + s] = std::cos(M_PI * angleY * straightAngleRec);
				element.ecf.angle[(6 * gpindex + 2) * SIMD::size + s] = std::cos(M_PI * angleZ * straightAngleRec);
				element.ecf.angle[(6 * gpindex + 3) * SIMD::size + s] = std::sin(M_PI * angleX * straightAngleRec);
				element.ecf.angle[(6 * gpindex + 4) * SIMD::size + s] = std::sin(M_PI * angleY * straightAngleRec);
				element.ecf.angle[(6 * gpindex + 5) * SIMD::size + s] = std::sin(M_PI * angleZ * straightAngleRec);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityRotationCylindric<nodes, gps, 2, edim, Physics>: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double cooX =    element.gpcoords[2 * gpindex + 0];
			double cooY =    element.gpcoords[2 * gpindex + 1];
			double centerX = element.ecf.angle [2 * gpindex + 0];
			double centerY = element.ecf.angle [2 * gpindex + 1];
			double rot = std::atan2(cooY - centerY, cooX - centerX);
			element.ecf.angle[2 * gpindex + 0] = std::cos(rot);
			element.ecf.angle[2 * gpindex + 1] = std::sin(rot);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				double cooX =    element.gpcoords[(2 * gpindex + 0) * SIMD::size + s];
				double cooY =    element.gpcoords[(2 * gpindex + 1) * SIMD::size + s];
				double centerX = element.ecf.angle [(2 * gpindex + 0) * SIMD::size + s];
				double centerY = element.ecf.angle [(2 * gpindex + 1) * SIMD::size + s];
				double rot = std::atan2(cooY - centerY, cooX - centerX);
				element.ecf.angle[(2 * gpindex + 0) * SIMD::size + s] = std::cos(rot);
				element.ecf.angle[(2 * gpindex + 1) * SIMD::size + s] = std::sin(rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityRotationCylindric<nodes, gps, 3, edim, Physics>: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double cooX =    element.gpcoords[3 * gpindex + 0];
			double cooY =    element.gpcoords[3 * gpindex + 1];
			double centerX = element.ecf.angle [6 * gpindex + 0];
			double centerY = element.ecf.angle [6 * gpindex + 1];
			double rot = std::atan2(cooY - centerY, cooX - centerX);
			element.ecf.angle[6 * gpindex + 0] = 1;
			element.ecf.angle[6 * gpindex + 1] = 1;
			element.ecf.angle[6 * gpindex + 2] = std::cos(rot);
			element.ecf.angle[6 * gpindex + 3] = 0;
			element.ecf.angle[6 * gpindex + 4] = 0;
			element.ecf.angle[6 * gpindex + 5] = std::sin(rot);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				double cooX =    element.gpcoords[(3 * gpindex + 0) * SIMD::size + s];
				double cooY =    element.gpcoords[(3 * gpindex + 1) * SIMD::size + s];
				double centerX = element.ecf.angle [(6 * gpindex + 0) * SIMD::size + s];
				double centerY = element.ecf.angle [(6 * gpindex + 1) * SIMD::size + s];
				double rot = std::atan2(cooY - centerY, cooX - centerX);
				element.ecf.angle[(6 * gpindex + 0) * SIMD::size + s] = 1;
				element.ecf.angle[(6 * gpindex + 1) * SIMD::size + s] = 1;
				element.ecf.angle[(6 * gpindex + 2) * SIMD::size + s] = std::cos(rot);
				element.ecf.angle[(6 * gpindex + 3) * SIMD::size + s] = 0;
				element.ecf.angle[(6 * gpindex + 4) * SIMD::size + s] = 0;
				element.ecf.angle[(6 * gpindex + 5) * SIMD::size + s] = std::sin(rot);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityRotationSpherical<nodes, gps, 3, edim, Physics>: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double x = element.gpcoords[3 * gpindex + 0] - element.ecf.angle [6 * gpindex + 0];
			double y = element.gpcoords[3 * gpindex + 1] - element.ecf.angle [6 * gpindex + 1];
			double z = element.gpcoords[3 * gpindex + 2] - element.ecf.angle [6 * gpindex + 2];
			double azimut = std::atan2(y, x);
			double r = std::sqrt(x * x + y * y + z * z);
			double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z * z + x * x), y);
			element.ecf.angle[6 * gpindex + 0] = 1;
			element.ecf.angle[6 * gpindex + 1] = std::cos(elevation);
			element.ecf.angle[6 * gpindex + 2] = std::cos(azimut);
			element.ecf.angle[6 * gpindex + 3] = 0;
			element.ecf.angle[6 * gpindex + 4] = std::sin(elevation);
			element.ecf.angle[6 * gpindex + 5] = std::sin(azimut);
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				double x = element.gpcoords[(3 * gpindex + 0) * SIMD::size + s] - element.ecf.angle [(6 * gpindex + 0) * SIMD::size + s];
				double y = element.gpcoords[(3 * gpindex + 1) * SIMD::size + s] - element.ecf.angle [(6 * gpindex + 1) * SIMD::size + s];
				double z = element.gpcoords[(3 * gpindex + 2) * SIMD::size + s] - element.ecf.angle [(6 * gpindex + 2) * SIMD::size + s];
				double azimut = std::atan2(y, x);
				double r = std::sqrt(x * x + y * y + z * z);
				double elevation = r < 1e-15 ? 0 : std::atan2(std::sqrt(z * z + x * x), y);
				element.ecf.angle[(6 * gpindex + 0) * SIMD::size + s] = 1;
				element.ecf.angle[(6 * gpindex + 1) * SIMD::size + s] = std::cos(elevation);
				element.ecf.angle[(6 * gpindex + 2) * SIMD::size + s] = std::cos(azimut);
				element.ecf.angle[(6 * gpindex + 3) * SIMD::size + s] = 0;
				element.ecf.angle[(6 * gpindex + 4) * SIMD::size + s] = std::sin(elevation);
				element.ecf.angle[(6 * gpindex + 5) * SIMD::size + s] = std::sin(azimut);
			}
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityRotationApply<nodes, gps, 2, edim, Physics>: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	// |cos , -sin| |c[0] , c[2]| | cos , sin|
	// |sin ,  cos| |c[3] , c[1]| |-sin , cos|
	// |cos * c[0] - sin * c[3] , cos * c[2] - sin * c[1]|  | cos , sin|
	// |sin * c[0] + cos * c[3] , sin * c[2] + cos * c[1]|  |-sin , cos|
	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double origin[4] = {
					element.ecf.conductivity[4 * gpindex + 0], element.ecf.conductivity[4 * gpindex + 1],
					element.ecf.conductivity[4 * gpindex + 2], element.ecf.conductivity[4 * gpindex + 3] };
			double cos = element.ecf.angle[2 * gpindex + 0];
			double sin = element.ecf.angle[2 * gpindex + 1];
			element.conductivity[4 * gpindex + 0] = (cos * origin[0] - sin * origin[2]) * cos - (cos * origin[1] - sin * origin[3]) * sin;
			element.conductivity[4 * gpindex + 1] = (cos * origin[0] - sin * origin[2]) * sin + (cos * origin[1] - sin * origin[3]) * cos;
			element.conductivity[4 * gpindex + 2] = (sin * origin[0] + cos * origin[2]) * cos - (sin * origin[1] + cos * origin[3]) * sin;
			element.conductivity[4 * gpindex + 3] = (sin * origin[0] + cos * origin[2]) * sin + (sin * origin[1] + cos * origin[3]) * cos;
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			SIMD origin0 = load(element.ecf.conductivity + (4 * gpindex + 0) * SIMD::size);
			SIMD origin1 = load(element.ecf.conductivity + (4 * gpindex + 1) * SIMD::size);
			SIMD origin2 = load(element.ecf.conductivity + (4 * gpindex + 2) * SIMD::size);
			SIMD origin3 = load(element.ecf.conductivity + (4 * gpindex + 3) * SIMD::size);
			SIMD cos = load(element.ecf.angle + (4 * gpindex + 0) * SIMD::size);
			SIMD sin = load(element.ecf.angle + (4 * gpindex + 1) * SIMD::size);

			SIMD res0 = (cos * origin0 - sin * origin2) * cos - (cos * origin1 - sin * origin3) * sin;
			SIMD res1 = (cos * origin0 - sin * origin2) * sin + (cos * origin1 - sin * origin3) * cos;
			SIMD res2 = (sin * origin0 + cos * origin2) * cos - (sin * origin1 + cos * origin3) * sin;
			SIMD res3 = (sin * origin0 + cos * origin2) * sin + (sin * origin1 + cos * origin3) * cos;

			store(element.conductivity + (4 * gpindex + 0) * SIMD::size, res0);
			store(element.conductivity + (4 * gpindex + 1) * SIMD::size, res1);
			store(element.conductivity + (4 * gpindex + 2) * SIMD::size, res2);
			store(element.conductivity + (4 * gpindex + 3) * SIMD::size, res3);
		}
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct ConductivityRotationApply<nodes, gps, 3, edim, Physics>: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double cos[3] { element.ecf.angle[6 * gpindex + 0], element.ecf.angle[6 * gpindex + 1], element.ecf.angle[6 * gpindex + 2] };
			double sin[3] { element.ecf.angle[6 * gpindex + 3], element.ecf.angle[6 * gpindex + 4], element.ecf.angle[6 * gpindex + 5] };
			double t[3][3] {
				{ cos[1] * cos[2]                           , cos[1] * sin[2]                           ,         -sin[1] },
				{ cos[2] * sin[0] * sin[1] - cos[0] * sin[2], cos[0] * cos[2] + sin[0] * sin[1] * sin[2], cos[1] * sin[0] },
				{ sin[0] * sin[2] + cos[0] * cos[2] * sin[1], cos[0] * sin[1] * sin[2] - cos[2] * sin[0], cos[0] * cos[1] }
			};

			double origin[9] = {
					element.ecf.conductivity[9 * gpindex + 0], element.ecf.conductivity[9 * gpindex + 1], element.ecf.conductivity[9 * gpindex + 2],
					element.ecf.conductivity[9 * gpindex + 3], element.ecf.conductivity[9 * gpindex + 4], element.ecf.conductivity[9 * gpindex + 5],
					element.ecf.conductivity[9 * gpindex + 6], element.ecf.conductivity[9 * gpindex + 7], element.ecf.conductivity[9 * gpindex + 8] };

			for (size_t i = 0; i < 3; ++i) {
				for (size_t j = 0; j < 3; ++j) {
					double _a = t[0][i] * origin[0] + t[1][i] * origin[3] + t[2][i] * origin[6];
					double _b = t[0][i] * origin[1] + t[1][i] * origin[4] + t[2][i] * origin[7];
					double _c = t[0][i] * origin[2] + t[1][i] * origin[5] + t[2][i] * origin[8];
					element.conductivity[9 * gpindex + 3 * i + j] = _a * t[0][j] + _b * t[1][j] + _c * t[2][j];
				}
			}
		}
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			SIMD cos0 = load(element.ecf.angle + 6 * gpindex + 0);
			SIMD cos1 = load(element.ecf.angle + 6 * gpindex + 1);
			SIMD cos2 = load(element.ecf.angle + 6 * gpindex + 2);
			SIMD sin0 = load(element.ecf.angle + 6 * gpindex + 3);
			SIMD sin1 = load(element.ecf.angle + 6 * gpindex + 4);
			SIMD sin2 = load(element.ecf.angle + 6 * gpindex + 5);

			SIMD t00 = cos1 * cos2;
			SIMD t01 = cos1 * sin2;
			SIMD t02 = -sin1;
			SIMD t10 = cos2 * sin0 * sin1 - cos0 * sin2;
			SIMD t11 = cos0 * cos2 + sin0 * sin1 * sin2;
			SIMD t12 = cos1 * sin0;
			SIMD t20 = sin0 * sin2 + cos0 * cos2 * sin1;
			SIMD t21 = cos0 * sin1 * sin2 - cos2 * sin0;
			SIMD t22 = cos0 * cos1;

			SIMD origin0 = load(element.ecf.conductivity + (9 * gpindex + 9) * SIMD::size);
			SIMD origin1 = load(element.ecf.conductivity + (9 * gpindex + 1) * SIMD::size);
			SIMD origin2 = load(element.ecf.conductivity + (9 * gpindex + 2) * SIMD::size);
			SIMD origin3 = load(element.ecf.conductivity + (9 * gpindex + 3) * SIMD::size);
			SIMD origin4 = load(element.ecf.conductivity + (9 * gpindex + 4) * SIMD::size);
			SIMD origin5 = load(element.ecf.conductivity + (9 * gpindex + 5) * SIMD::size);
			SIMD origin6 = load(element.ecf.conductivity + (9 * gpindex + 6) * SIMD::size);
			SIMD origin7 = load(element.ecf.conductivity + (9 * gpindex + 7) * SIMD::size);
			SIMD origin8 = load(element.ecf.conductivity + (9 * gpindex + 8) * SIMD::size);

			SIMD _a = t00 * origin0 + t10 * origin3 + t20 * origin6;
			SIMD _b = t00 * origin1 + t10 * origin4 + t20 * origin7;
			SIMD _c = t00 * origin2 + t10 * origin5 + t20 * origin8;

			SIMD res = _a * t00 + _b * t10 + _c * t20;
			store(element.conductivity + (9 * gpindex + 0) * SIMD::size, res);

			res = _a * t01 + _b * t11 + _c * t21;
			store(element.conductivity + (9 * gpindex + 1) * SIMD::size, res);

			res = _a * t02 + _b * t12 + _c * t22;
			store(element.conductivity + (9 * gpindex + 2) * SIMD::size, res);

			_a = t01 * origin0 + t11 * origin3 + t21 * origin6;
			_b = t01 * origin1 + t11 * origin4 + t21 * origin7;
			_c = t01 * origin2 + t11 * origin5 + t21 * origin8;

			res = _a * t00 + _b * t10 + _c * t20;
			store(element.conductivity + (9 * gpindex + 3) * SIMD::size, res);

			res = _a * t01 + _b * t11 + _c * t21;
			store(element.conductivity + (9 * gpindex + 4) * SIMD::size, res);

			res = _a * t02 + _b * t12 + _c * t22;
			store(element.conductivity + (9 * gpindex + 5) * SIMD::size, res);

			_a = t02 * origin0 + t12 * origin3 + t22 * origin6;
			_b = t02 * origin1 + t12 * origin4 + t22 * origin7;
			_c = t02 * origin2 + t12 * origin5 + t22 * origin8;

			res = _a * t00 + _b * t10 + _c * t20;
			store(element.conductivity + (9 * gpindex + 6) * SIMD::size, res);

			res = _a * t01 + _b * t11 + _c * t21;
			store(element.conductivity + (9 * gpindex + 7) * SIMD::size, res);

			res = _a * t02 + _b * t12 + _c * t22;
			store(element.conductivity + (9 * gpindex + 8) * SIMD::size, res);
		}
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct ConductivityRotationSkip: ConductivityRotation, Physics {
	using ConductivityRotation::ConductivityRotation;

	void sisd(typename Physics::Element &element)
	{
		memcpy(element.conductivity, element.ecf.conductivity, ndim * ndim * gps * sizeof(double));
	}

	void simd(typename Physics::Element &element)
	{
		memcpy(element.conductivity, element.ecf.conductivity, ndim * ndim * gps * sizeof(double) * SIMD::size);
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_CONDUCTIVITY_H_ */
