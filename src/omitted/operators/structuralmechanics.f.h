
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_F_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_F_H_

#include <analysis/assembler/subkernel/operator.h>

namespace espreso {

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct Acceleration;

template <size_t nodes, size_t gps, class Physics>
struct Acceleration<nodes, gps, 2, 2, StructuralMechanicsElementType::SYMMETRIC_PLANE, Physics>: ActionOperator, Physics {
	const char* name() const { return "Acceleration"; }

	OutputParameterIterator rhs;

	Acceleration(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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
	const char* name() const { return "Acceleration"; }

	OutputParameterIterator rhs;

	Acceleration(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = rhs.data;
		for (size_t n = 0; n < nodes; ++n) {
//			SIMD fx = load(out + (0 * nodes + n) * SIMD::size);
			SIMD fy = load(out + (1 * nodes + n) * SIMD::size);
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * element.ecf.density[gp] * load1(element.N[gp][n]);
//				fx = fx + scale * element.ecf.acceleration[gp][0];
				fy = fy + scale * element.ecf.acceleration[gp][1];
			}
//			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct Acceleration<nodes, gps, 3, 3, StructuralMechanicsElementType::SYMMETRIC_VOLUME, Physics>: ActionOperator, Physics {
	const char* name() const { return "Acceleration"; }

	OutputParameterIterator rhs;

	Acceleration(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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
	const char* name() const { return "AngularVelocity"; }

	OutputParameterIterator rhs;

	AngularVelocity(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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
	const char* name() const { return "AngularVelocity"; }

	OutputParameterIterator rhs;

	AngularVelocity(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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
			for (size_t gp = 0; gp < gps; ++gp) {
				SIMD scale = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * element.ecf.density[gp] * load1(element.N[gp][n]) * element.ecf.angularVelocity[gp] * element.ecf.angularVelocity[gp];
				fx = fx + scale * element.gpcoords[gp][0];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct AngularVelocity<nodes, gps, 3, 3, StructuralMechanicsElementType::SYMMETRIC_VOLUME, Physics>: ActionOperator, Physics {
	const char* name() const { return "AngularVelocity"; }

	OutputParameterIterator rhs;

	AngularVelocity(size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct NormalPressure;

template <size_t nodes, size_t gps, class Physics>
struct NormalPressure<nodes, gps, 2, 1, StructuralMechanicsElementType::EDGE, Physics>: ActionOperator, Physics {
	const char* name() const { return "NormalPressure"; }

	OutputParameterIterator rhs;

	NormalPressure(size_t region, size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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
				SIMD pressure = element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.normalPressure[gp];
				fx = fx + pressure * element.normal[gp][0];
				fy = fy + pressure * element.normal[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct NormalPressure<nodes, gps, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC, Physics>: ActionOperator, Physics {
	const char* name() const { return "NormalPressure"; }

	OutputParameterIterator rhs;

	NormalPressure(size_t region, size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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
				SIMD pressure = element.det[gp] * load1(element.w[gp]) * load1(2 * M_PI) * element.gpcoords[gp][0] * load1(element.N[gp][n]) * element.ecf.normalPressure[gp];
				fx = fx + pressure * element.normal[gp][0];
				fy = fy + pressure * element.normal[gp][1];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct NormalPressure<nodes, gps, 3, 2, StructuralMechanicsElementType::FACE, Physics>: ActionOperator, Physics {
	const char* name() const { return "NormalPressure"; }

	OutputParameterIterator rhs;

	NormalPressure(size_t region, size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
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
				SIMD pressure = element.det[gp] * load1(element.w[gp]) * load1(element.N[gp][n]) * element.ecf.normalPressure[gp];
				fx = fx + pressure * element.normal[gp][0];
				fy = fy + pressure * element.normal[gp][1];
				fz = fz + pressure * element.normal[gp][2];
			}
			store(out + (0 * nodes + n) * SIMD::size, fx);
			store(out + (1 * nodes + n) * SIMD::size, fy);
			store(out + (2 * nodes + n) * SIMD::size, fz);
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, class Physics>
struct NormalPressure<nodes, gps, 3, 1, StructuralMechanicsElementType::EDGE, Physics>: ActionOperator, Physics {
	const char* name() const { return "NormalPressure"; }

	OutputParameterIterator rhs;

	NormalPressure(size_t region, size_t interval, ParameterData &rhs)
	: rhs(rhs, interval)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void move(int n)
	{
		rhs += n;
	}

	void simd(typename Physics::Element &element)
	{

	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_F_H_ */
