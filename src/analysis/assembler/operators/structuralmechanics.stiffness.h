
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_STIFFNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_STIFFNESS_H_

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

}



#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_STRUCTURALMECHANICS_STIFFNESS_H_ */
