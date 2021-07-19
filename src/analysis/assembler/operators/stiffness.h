
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

struct Stiffness: public ActionOperator {
	Stiffness(
			int interval,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &conductivity,
			const ParameterData &xi,
			const ParameterData &thickness,
			ParameterData &stiffness)
	: ActionOperator(interval, stiffness.isconst[interval], stiffness.update[interval]),
	  dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  conductivity(conductivity, interval),
	  xi(xi, interval),
	  thickness(thickness, interval),
	  stiffness(stiffness, interval)
	{
		if (update) {
			std::fill((stiffness.data->begin() + interval)->data(), (stiffness.data->begin() + interval + 1)->data(), 0);
		}
	}

	InputParameterIterator dND, weight, determinant, conductivity, xi, thickness;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant; ++conductivity; ++xi; ++thickness;
		++stiffness;
	}

	Stiffness& operator+=(const size_t rhs)
	{
		dND += rhs; determinant += rhs; conductivity += rhs; xi += rhs; thickness += rhs;
		stiffness += rhs;
		return *this;
	}

	void reset()
	{

	}
};

template<size_t nodes, size_t gps>
struct Stiffness2DHeatIsotropic: public Stiffness {
	using Stiffness::Stiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN2M2N<nodes>(thickness[gpindex] * xi[gpindex] * determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness2DHeat: public Stiffness {
	using Stiffness::Stiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN2M22M2N<nodes>(thickness[gpindex] * xi[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 4 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness3DHeatIsotropic: public Stiffness {
	using Stiffness::Stiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN3M3N<nodes>(xi[gpindex] * determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness3DHeat: public Stiffness {
	using Stiffness::Stiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN3M33M3N<nodes>(xi[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_STIFFNESS_H_ */
