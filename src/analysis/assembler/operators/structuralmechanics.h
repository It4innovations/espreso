
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

struct StructuralMechanicsStiffness: public ActionOperator {
	StructuralMechanicsStiffness(
			int interval,
			const ParameterData &N,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &coordinates,
			const ParameterData &elasticity,
			const ParameterData &thickness,
			ParameterData &stiffness)
	: N(N, interval),
	  dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  coordinates(coordinates, interval),
	  elasticity(elasticity, interval),
	  thickness(thickness, interval),
	  stiffness(stiffness, interval)
	{

	}

	InputParameterIterator N, dND, weight, determinant, coordinates, elasticity, thickness;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant; ++coordinates; ++elasticity; ++thickness;
		++stiffness;
	}

	void move(int n)
	{
		dND += n; determinant += n; coordinates += n; elasticity += n; thickness += n;
		stiffness += n;
	}
};

template<size_t nodes, size_t gps>
struct StiffnessPlane: public StructuralMechanicsStiffness {
	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDDXDY33DXDYN<nodes>(determinant[gpindex] * weight[gpindex], elasticity.data + 9 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct StiffnessPlaneWithThickness: public StructuralMechanicsStiffness {
	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDDXDY33DXDYN<nodes>(determinant[gpindex] * weight[gpindex] * thickness[gpindex], elasticity.data + 9 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct StiffnessAxisymmetric: public StructuralMechanicsStiffness {
	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			double coo[nodes];
			for (int n = 0; n < nodes; ++n) {
				coo[n] = N[gpindex * nodes + n] / coordinates[2 * gpindex];
			}
			ADDDXDYCOO44DXDYN<nodes>(determinant[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex], elasticity.data + 16 * gpindex, dND.data + 2 * nodes * gpindex, coo, stiffness.data);
		}
	}
};


template<size_t nodes, size_t gps>
struct Stiffness3DElasticity: public StructuralMechanicsStiffness {
	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDDXDYDZ66DXDYDZN<nodes>(determinant[gpindex] * weight[gpindex], elasticity.data + 36 * gpindex, dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

struct Acceleration: public ActionOperator {
	Acceleration(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &density, const ParameterData &coordinates, const ParameterData &thickness, const ParameterData &acceleration, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  density(density, interval),
	  coordinates(coordinates, interval),
	  thickness(thickness, interval),
	  acceleration(acceleration, interval),
	  rhs(rhs, interval)
	{

	}

	InputParameterIterator N, weight, J, density, coordinates, thickness;
	InputParameterIterator acceleration;
	OutputParameterIterator rhs;

	void operator++()
	{
		++weight; ++J; ++density; ++thickness;
		++acceleration;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n; density += n; thickness += n;
		acceleration += n;
		rhs += n;
	}
};

template<size_t nodes, size_t gps>
struct AccelerationPlane: public Acceleration {
	using Acceleration::Acceleration;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[2 * gpindex + 0];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[2 * gpindex + 1];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct AccelerationPlaneWithThickness: public Acceleration {
	using Acceleration::Acceleration;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * thickness[gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[2 * gpindex + 0];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * thickness[gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[2 * gpindex + 1];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct AccelerationAxisymmetric: public Acceleration {
	using Acceleration::Acceleration;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[2 * gpindex + 0];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[2 * gpindex + 1];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct Acceleration3D: public Acceleration {
	using Acceleration::Acceleration;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[0 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[3 * gpindex + 0];
				rhs[1 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[3 * gpindex + 1];
				rhs[2 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * acceleration[3 * gpindex + 2];
			}
		}
	}
};

struct AngularVelocity: public ActionOperator {
	AngularVelocity(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &density, const ParameterData &thickness, const ParameterData &coordinates, const ParameterData &angularVelocity, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  density(density, interval),
	  thickness(thickness, interval),
	  coordinates(coordinates, interval),
	  angularVelocity(angularVelocity, interval),
	  rhs(rhs, interval)
	{

	}

	InputParameterIterator N, weight, J, density, thickness;
	InputParameterIterator coordinates, angularVelocity;
	OutputParameterIterator rhs;

	void operator++()
	{
		++weight; ++J; ++density; ++thickness;
		++coordinates; ++angularVelocity;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n; density += n; thickness += n;
		coordinates += n; angularVelocity += n;
		rhs += n;
	}
};

template<size_t nodes, size_t gps>
struct AngularVelocityPlane: public AngularVelocity {
	using AngularVelocity::AngularVelocity;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[2 * gpindex + 0] * angularVelocity[3 * gpindex + 2] * angularVelocity[3 * gpindex + 2];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[2 * gpindex + 1] * angularVelocity[3 * gpindex + 2] * angularVelocity[3 * gpindex + 2];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct AngularVelocityPlaneWithThickness: public AngularVelocity {
	using AngularVelocity::AngularVelocity;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * thickness[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[2 * gpindex + 0] * angularVelocity[3 * gpindex + 2] * angularVelocity[3 * gpindex + 2];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * thickness[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[2 * gpindex + 1] * angularVelocity[3 * gpindex + 2] * angularVelocity[3 * gpindex + 2];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct AngularVelocityAxisymmetric: public AngularVelocity {
	using AngularVelocity::AngularVelocity;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[2 * gpindex + 0] * angularVelocity[3 * gpindex + 1] * angularVelocity[3 * gpindex + 1];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[2 * gpindex + 1] * angularVelocity[3 * gpindex + 1] * angularVelocity[3 * gpindex + 1];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct AngularVelocity3D: public AngularVelocity {
	using AngularVelocity::AngularVelocity;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				// angular velocity X
				rhs[1 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[3 * gpindex + 1] * angularVelocity[3 * gpindex + 0] * angularVelocity[3 * gpindex + 0];
				rhs[2 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[3 * gpindex + 2] * angularVelocity[3 * gpindex + 0] * angularVelocity[3 * gpindex + 0];

				// angular velocity Y
				rhs[0 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[3 * gpindex + 0] * angularVelocity[3 * gpindex + 1] * angularVelocity[3 * gpindex + 1];
				rhs[2 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[3 * gpindex + 2] * angularVelocity[3 * gpindex + 1] * angularVelocity[3 * gpindex + 1];

				// angular velocity Z
				rhs[0 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[3 * gpindex + 0] * angularVelocity[3 * gpindex + 2] * angularVelocity[3 * gpindex + 2];
				rhs[1 * nodes + n] += J[gpindex] * weight[gpindex] * density[gpindex] * N[gpindex * nodes + n] * coordinates[3 * gpindex + 1] * angularVelocity[3 * gpindex + 2] * angularVelocity[3 * gpindex + 2];
			}
		}
	}
};

struct NormalPressure: public ActionOperator {
	NormalPressure(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &normal, const ParameterData &coordinates, const ParameterData &thickness, const ParameterData &normalPressure, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  normal(normal, interval),
	  coordinates(coordinates, interval),
//	  thickness(thickness, interval), // TODO
	  normalPressure(normalPressure, interval),
	  rhs(rhs, interval)
	{ }

	InputParameterIterator N, weight, J, normal, coordinates;
//	InputParameterIterator thickness;
	InputParameterIterator normalPressure;
	OutputParameterIterator rhs;

	void operator++()
	{
		++weight; ++J, ++normal; ++coordinates;
//		++thickness;
		++normalPressure;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n; normal += n; coordinates += n;
//		thickness += n;
		normalPressure += n;
		rhs += n;
	}
};

template<size_t nodes, size_t gps>
struct NormalPressurePlane: public NormalPressure {
	using NormalPressure::NormalPressure;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct NormalPressurePlaneWithThickness: public NormalPressure {
	using NormalPressure::NormalPressure;

	void operator()()
	{
		// TODO: thickness
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * normal[2 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct NormalPressureAxisymmetric: public NormalPressure {
	using NormalPressure::NormalPressure;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[        n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * normal[2 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
				rhs[nodes + n] += J[gpindex] * weight[gpindex] * 2 * M_PI * coordinates[2 * gpindex] * normal[2 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct NormalPressure3D: public NormalPressure {
	using NormalPressure::NormalPressure;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[0 * nodes + n] += J[gpindex] * weight[gpindex] * normal[3 * gpindex + 0] * normalPressure[gpindex] * N[gpindex * nodes + n];
				rhs[1 * nodes + n] += J[gpindex] * weight[gpindex] * normal[3 * gpindex + 1] * normalPressure[gpindex] * N[gpindex * nodes + n];
				rhs[2 * nodes + n] += J[gpindex] * weight[gpindex] * normal[3 * gpindex + 2] * normalPressure[gpindex] * N[gpindex * nodes + n];
			}
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_H_ */
