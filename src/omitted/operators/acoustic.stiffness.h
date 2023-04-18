
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_STIFFNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_STIFFNESS_H_

#include <analysis/assembler/subkernel/operator.h>
#include "analysis/assembler/parameter.h"

#include <iostream>
#include <complex>

namespace espreso {

struct AcousticStiffness: public ActionOperator {
	const char* name() const { return "AcousticStiffness"; }

	AcousticStiffness(
			int interval,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &density,
			ParameterData &stiffness)
	: dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  density(density, interval),
	  stiffness(stiffness, interval)
	{

	}

	InputParameterIterator dND, weight, determinant, density;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant;;
		++stiffness;
		++density;
	}

	void move(int n)
	{
		dND += n; determinant += n;
		stiffness += n;
		density += n;
	}

	AcousticStiffness& operator+=(const size_t rhs)
	{
		dND += rhs; determinant += rhs;;
		stiffness += rhs;
		return *this;
	}
};


template<size_t nodes, size_t gps>
struct Stiffness2DAcoustic: public AcousticStiffness {
	using AcousticStiffness::AcousticStiffness;

	void operator()()
	{
		std::fill(stiffness.data, stiffness.data + stiffness.inc, 0);
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDMN2M2N<nodes>(determinant[gpindex] * weight[gpindex] / density[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness3DAcoustic: public AcousticStiffness {
	using AcousticStiffness::AcousticStiffness;

	void operator()()
	{
		std::fill(stiffness.data, stiffness.data + stiffness.inc, 0);
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDMN3M3N<nodes>(determinant[gpindex] * weight[gpindex] / density[gpindex], dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

template <size_t nodes, size_t gps>
struct AcousticMass: public ActionOperator {
	const char* name() const { return "AcousticMass"; }

	InputParameterIterator N, weight, determinant, density, speed_of_sound;
	OutputParameterIterator mass;

	AcousticMass(
		int interval,
		const ParameterData &N,
		const ParameterData &weight,
		const ParameterData &determinant,
		const ParameterData &density,
		const ParameterData &speed_of_sound,
		ParameterData &mass)
	: N(N, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  density(density, interval),
	  speed_of_sound(speed_of_sound, interval),
	  mass(mass, interval)
	{

	}

	void operator++()
	{
		++N; ++determinant;
		++mass;
		++density;
	}

	void move(int n)
	{
		N += n; determinant += n;
		mass += n;
		density += n;
	}

	void operator()()
	{
		std::fill(mass.data, mass.data + mass.inc, 0);
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDMN1M1N<nodes>(determinant[gpindex] * weight[gpindex] / (density[gpindex] * speed_of_sound[gpindex] * speed_of_sound[gpindex]), N.data + nodes * gpindex, mass.data);
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_STIFFNESS_H_ */
