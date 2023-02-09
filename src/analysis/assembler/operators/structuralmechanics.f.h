
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_F_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_F_H_

#include "analysis/assembler/operator.h"

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct Acceleration;

template <size_t nodes, size_t gps, class Physics>
struct Acceleration<nodes, gps, 2, 2, StructuralMechanicsElementType::SYMMETRIC_PLANE, Physics>: ActionOperator, Physics {
	OutputParameterIterator rhs;

	Acceleration(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]);
				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct Acceleration<nodes, gps, 2, 2, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC, Physics>: ActionOperator, Physics {
	OutputParameterIterator rhs;

	Acceleration(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * element.ecf.density[gp] * load1(element.N[gp][n]);
				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct Acceleration<nodes, gps, 3, 3, StructuralMechanicsElementType::SYMMETRIC_VOLUME, Physics>: ActionOperator, Physics {
	OutputParameterIterator rhs;

	Acceleration(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			SIMD fz = load(out + (2 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]);
				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
				fz = fz + scale * element.ecf.acceleration[gp][2];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
			store(out + (2 * nodes + n) * SIMD::size, fz);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct AngularVelocity;

template <size_t nodes, size_t gps, class Physics>
struct AngularVelocity<nodes, gps, 2, 2, StructuralMechanicsElementType::SYMMETRIC_PLANE, Physics>: ActionOperator, Physics {
	OutputParameterIterator rhs;

	AngularVelocity(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.ecf.thickness[gp] * element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]) * element.ecf.angularVelocity[gp] * element.ecf.angularVelocity[gp];
				fx = fx + scale * element.gpcoords[gp][0];
				fy = fy + scale * element.gpcoords[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct AngularVelocity<nodes, gps, 2, 2, StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC, Physics>: ActionOperator, Physics {
	OutputParameterIterator rhs;

	AngularVelocity(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * element.ecf.density[gp] * load1(element.N[gp][n]) * element.ecf.angularVelocity[gp] * element.ecf.angularVelocity[gp];
				fx = fx + scale * element.gpcoords[gp][0];
				fy = fy + scale * element.gpcoords[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct AngularVelocity<nodes, gps, 3, 3, StructuralMechanicsElementType::SYMMETRIC_VOLUME, Physics>: ActionOperator, Physics {
	OutputParameterIterator rhs;

	AngularVelocity(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			SIMD fz = load(out + (2 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * element.ecf.density[gp] * load1(element.N[gp][n]);
				SIMD vx = element.ecf.angularVelocity[gp][0] * element.ecf.angularVelocity[gp][0];
				SIMD vy = element.ecf.angularVelocity[gp][1] * element.ecf.angularVelocity[gp][1];
				SIMD vz = element.ecf.angularVelocity[gp][2] * element.ecf.angularVelocity[gp][2];
				fx = fx + scale * element.gpcoords[gp][0] * (vy + vz);
				fy = fy + scale * element.gpcoords[gp][1] * (vx + vz);
				fz = fz + scale * element.gpcoords[gp][2] * (vx + vy);
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
			store(out + (2 * nodes + n) * SIMD::size, fz);
		}
		move(SIMD::size);
	}
};

//struct NormalPressure: public ActionOperator {
//	NormalPressure(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &normal, const ParameterData &coordinates, const ParameterData &thickness, const ParameterData &normalPressure, ParameterData &rhs)
//	: N(N, interval),
//	  weight(weight, interval),
//	  J(J, interval),
//	  normal(normal, interval),
//	  coordinates(coordinates, interval),
////	  thickness(thickness, interval), // TODO
//	  normalPressure(normalPressure, interval),
//	  rhs(rhs, interval)
//	{ }
//
//	InputParameterIterator N, weight, J, normal, coordinates;
////	InputParameterIterator thickness;
//	InputParameterIterator normalPressure;
//	OutputParameterIterator rhs;
//
//	void operator++()
//	{
//		++weight; ++J, ++normal; ++coordinates;
////		++thickness;
//		++normalPressure;
//		++rhs;
//	}
//
//	void move(int n)
//	{
//		weight += n; J += n; normal += n; coordinates += n;
////		thickness += n;
//		normalPressure += n;
//		rhs += n;
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct NormalPressurePlane: public NormalPressure {
//	using NormalPressure::NormalPressure;
//
//	void operator()()
//	{
//		for (size_t n = 0; n < nodes; ++n) {
//			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//				rhs[        n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
//				rhs[nodes + n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
//			}
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct NormalPressurePlaneWithThickness: public NormalPressure {
//	using NormalPressure::NormalPressure;
//
//	void operator()()
//	{
//		// TODO: thickness
//		for (size_t n = 0; n < nodes; ++n) {
//			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//				rhs[        n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
//				rhs[nodes + n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
//			}
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct NormalPressureAxisymmetric: public NormalPressure {
//	using NormalPressure::NormalPressure;
//
//	void operator()()
//	{
//		for (size_t n = 0; n < nodes; ++n) {
//			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//				rhs[        n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * normal[2 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
//				rhs[nodes + n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * normal[2 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
//			}
//		}
//	}
//};
//
//template<size_t nodes, size_t gps>
//struct NormalPressure3D: public NormalPressure {
//	using NormalPressure::NormalPressure;
//
//	void operator()()
//	{
//		for (size_t n = 0; n < nodes; ++n) {
//			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//				rhs[0 * nodes + n] += J[gpindex] * weight[gpindex] * normal[3 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
//				rhs[1 * nodes + n] += J[gpindex] * weight[gpindex] * normal[3 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
//				rhs[2 * nodes + n] += J[gpindex] * weight[gpindex] * normal[3 * gpindex + 2] * normalPressure[gpindex] * N[gpindex * nodes + n];
//			}
//		}
//	}
//};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_F_H_ */
