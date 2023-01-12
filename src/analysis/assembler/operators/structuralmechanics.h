
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

struct StructuralMechanicsStiffness: public ActionOperator {
	StructuralMechanicsStiffness(
			int interval,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &elasticity,
			const ParameterData &thickness,
			ParameterData &stiffness)
	: dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  elasticity(elasticity, interval),
	  thickness(thickness, interval),
	  stiffness(stiffness, interval)
	{

	}

	InputParameterIterator dND, weight, determinant, elasticity,thickness;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant; ++elasticity; ++thickness;
		++stiffness;
	}

	void move(int n)
	{
		dND += n; determinant += n; elasticity += n; thickness += n;
		stiffness += n;
	}

	StructuralMechanicsStiffness& operator+=(const size_t rhs)
	{
		dND += rhs; determinant += rhs; elasticity += rhs; thickness += rhs;
		stiffness += rhs;
		return *this;
	}
};

template<size_t nodes, size_t gps>
struct Stiffness2DPlane: public StructuralMechanicsStiffness {
	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDDXDY33DXDYN<nodes>(determinant[gpindex] * weight[gpindex], elasticity.data + 9 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness2DPlaneWithThickness: public StructuralMechanicsStiffness {
	using StructuralMechanicsStiffness::StructuralMechanicsStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDDXDY33DXDYN<nodes>(determinant[gpindex] * weight[gpindex] * thickness[gpindex], elasticity.data + 9 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
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

template<size_t nodes, size_t gps>
struct StructuralMechanicsRHS: public ActionOperator {
	StructuralMechanicsRHS(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &heatSource, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  heatSource(heatSource, interval),
	  rhs(rhs, interval)
	{

	}

	InputParameterIterator N, weight, J;
	InputParameterIterator heatSource;
	OutputParameterIterator rhs;

	void operator()()
	{
		std::fill(rhs.data, rhs.data + rhs.inc, 0);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs.data[n] += J.data[gpindex] * weight.data[gpindex] * heatSource.data[gpindex] * N.data[gpindex * nodes + n];
			}
		}
	}

	void operator++()
	{
		++weight; ++J;
		++heatSource;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n;
		heatSource += n;
		rhs += n;
	}
};

struct NormalPressure: public ActionOperator {
	NormalPressure(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &normal, const ParameterData &thickness, const ParameterData &normalPressure, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  normal(normal, interval),
//	  thickness(thickness, interval), // TODO
	  normalPressure(normalPressure, interval),
	  rhs(rhs, interval)
	{ }

	InputParameterIterator N, weight, J, normal;
//	InputParameterIterator thickness;
	InputParameterIterator normalPressure;
	OutputParameterIterator rhs;

	void operator++()
	{
		++weight; ++J, ++normal;
//		++thickness;
		++normalPressure;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n; normal += n;
//		thickness += n;
		normalPressure += n;
		rhs += n;
	}
};

template<size_t nodes, size_t gps>
struct NormalPressure2D: public NormalPressure {
	using NormalPressure::NormalPressure;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs.data[        n] += J.data[gpindex] * weight.data[gpindex] * normal.data[2 * gpindex + 0] * normalPressure.data[gpindex] * N.data[gpindex * nodes + n];
				rhs.data[nodes + n] += J.data[gpindex] * weight.data[gpindex] * normal.data[2 * gpindex + 1] * normalPressure.data[gpindex] * N.data[gpindex * nodes + n];
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
				rhs.data[0 * nodes + n] += J.data[gpindex] * weight.data[gpindex] * normal.data[3 * gpindex + 0] * normalPressure.data[gpindex] * N.data[gpindex * nodes + n];
				rhs.data[1 * nodes + n] += J.data[gpindex] * weight.data[gpindex] * normal.data[3 * gpindex + 1] * normalPressure.data[gpindex] * N.data[gpindex * nodes + n];
				rhs.data[2 * nodes + n] += J.data[gpindex] * weight.data[gpindex] * normal.data[3 * gpindex + 2] * normalPressure.data[gpindex] * N.data[gpindex * nodes + n];
			}
		}
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_H_ */
