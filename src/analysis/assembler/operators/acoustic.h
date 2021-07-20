
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

struct AcousticStiffness: public ActionOperator {
	AcousticStiffness(
			int interval,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			ParameterData &stiffness)
	: ActionOperator(interval, stiffness.isconst[interval], stiffness.update[interval]),
	  dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  stiffness(stiffness, interval)
	{

	}

	InputParameterIterator dND, weight, determinant;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant;;
		++stiffness;
	}

	AcousticStiffness& operator+=(const size_t rhs)
	{
		dND += rhs; determinant += rhs;;
		stiffness += rhs;
		return *this;
	}

	void reset()
	{

	}
};


template<size_t nodes, size_t gps>
struct Stiffness2DAcoustic: public AcousticStiffness {
	using AcousticStiffness::AcousticStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN2M2N<nodes>(determinant[gpindex] * weight[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness3DAcoustic: public AcousticStiffness {
	using AcousticStiffness::AcousticStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN3M3N<nodes>(determinant[gpindex] * weight[gpindex], dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

template <size_t nodes, size_t gps>
struct AcousticMass: public ActionOperator {
	InputParameterIterator N, weight, determinant;
	OutputParameterIterator mass;

	AcousticMass(
		int interval,
		const ParameterData &N,
		const ParameterData &weight,
		const ParameterData &determinant,
		ParameterData &mass)
	: ActionOperator(interval, mass.isconst[interval], mass.update[interval]),
	  N(N, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  mass(mass, interval)
	{

	}

	void operator++()
	{
		++N; ++determinant;
		++mass;
	}

	void operator()()
	{
		// double kappa = omega / speedOfSound
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN1M1N<nodes>(determinant[gpindex] * weight[gpindex] /*  * kappa[gpindex] * kappa[gpindex] */, N.data + nodes * gpindex, mass.data);
		}
	}

	void reset()
	{

	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_H_ */
