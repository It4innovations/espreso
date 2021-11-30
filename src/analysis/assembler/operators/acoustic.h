
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

#include <iostream>

namespace espreso {

struct AcousticDipole: public ActionOperator {
	AcousticDipole(
			int interval,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &density,
			const ParameterData &q,
			ParameterData &dipole)

	: dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  density(density, interval),
	  q(q, interval),
	  dipole(dipole, interval)
	{
	}

	InputParameterIterator dND, weight, determinant, density, q;
	OutputParameterIterator dipole;

	void operator++()
	{
		++dND;
		++determinant;;
		++dipole;
		++density;
		++q;
	}

	void move(int n)
	{
		dND += n;
		determinant += n;
		dipole += n;
		density += n;
		q += n;
	}
};

template<size_t nodes, size_t gps>
struct AcousticDipole2D: public AcousticDipole {
	using AcousticDipole::AcousticDipole;

	void operator()()
	{
		std::fill(dipole.data, dipole.data + dipole.inc, 0);

		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				dipole.data[n] += (determinant[gpindex] * weight[gpindex] / density[gpindex]) * (q[2 * gpindex + 0] * dND[gpindex * nodes * 2 + 0 * nodes + n] + q[2 * gpindex + 1] * dND[gpindex * nodes * 2 + 1 * nodes + n]);
			}
		}
	}
};

template<size_t nodes, size_t gps>
struct AcousticDipole3D: public AcousticDipole {
	using AcousticDipole::AcousticDipole;
	
	void operator()()
	{
		std::fill(dipole.data, dipole.data + dipole.inc, 0);

		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				dipole.data[n] += (determinant[gpindex] * weight[gpindex] / density[gpindex]) * (q[3 * gpindex + 0] * dND[gpindex * nodes * 3 + 0 * nodes + n] + q[3 * gpindex + 1] * dND[gpindex * nodes * 3 + 1 * nodes + n] + q[3 * gpindex + 2] * dND[gpindex * nodes * 3 + 2 * nodes + n]);
			}
		}
	}
};


struct AcousticStiffness: public ActionOperator {
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
			ADDMN2M2N<nodes>(determinant[gpindex] * weight[gpindex] / density[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
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
			ADDMN3M3N<nodes>(determinant[gpindex] * weight[gpindex] / density[gpindex], dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

template <size_t nodes, size_t gps>
struct AcousticMass: public ActionOperator {
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
			ADDMN1M1N<nodes>(determinant[gpindex] * weight[gpindex] / (density[gpindex] * speed_of_sound[gpindex] * speed_of_sound[gpindex]), N.data + nodes * gpindex, mass.data);
		}
	}
};


template <size_t nodes, size_t gps>
struct AcousticQ: public ActionOperator {
	AcousticQ(int interval, double area, const ParameterData &g,  ParameterData &impedance, ParameterData &q)
	: area(area),
	  g(g, interval),
	  impedance(impedance, interval),
	  q(q, interval)
	{

	}

	void operator()()
	{
		std::fill(q.data, q.data + q.inc, 0);
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			q.data[gpindex] += g.data[gpindex];
		}
	}

	void operator++()
	{
		++g;
		++q;
	}

	void move(int n)
	{
		g += n;
		q += n;
	}

	double area;
	InputParameterIterator g, impedance;
	OutputParameterIterator q;
};

template <size_t nodes, size_t gps>
struct AcousticRHS2D: public ActionOperator {
	AcousticRHS2D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &q, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  q(q, interval),
	  rhs(rhs, interval)
	{ }

	InputParameterIterator N, weight, J;
	InputParameterIterator q;
	OutputParameterIterator rhs;

	void operator()()
	{
		std::fill(rhs.data, rhs.data + rhs.inc, 0);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs.data[n] += J.data[gpindex] * weight.data[gpindex] * q.data[gpindex] * N.data[gpindex * nodes + n];
			}
		}
	}

	void operator++()
	{
		++weight; ++J;
		++q;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n;
		q += n;
		rhs += n;
	}
};

template <size_t nodes, size_t gps>
struct AcousticRHS3D: public ActionOperator {
	AcousticRHS3D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &q, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  q(q, interval),
	  rhs(rhs, interval)
	{

	}

	InputParameterIterator N, weight, J;
	InputParameterIterator q;
	OutputParameterIterator rhs;

	void operator()()
	{
		std::fill(rhs.data, rhs.data + rhs.inc, 0);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs.data[n] += J.data[gpindex] * weight.data[gpindex] * q.data[gpindex] * N.data[gpindex * nodes + n];
			}
		}
	}

	void operator++()
	{
		++weight; ++J;
		++q;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n;
		q += n;
		rhs += n;
	}
};

template <size_t nodes, size_t gps>
struct AcousticsBoundaryMass: public ActionOperator {
	InputParameterIterator N, weight, determinant, impedance;
	OutputParameterIterator boundaryMass;

	AcousticsBoundaryMass(
		int interval,
		const ParameterData &N,
		const ParameterData &weight,
		const ParameterData &determinant,
		const ParameterData &impedance,
		ParameterData &boundaryMass)
	: N(N, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  impedance(impedance, interval),
	  boundaryMass(boundaryMass, interval)
	{

	}

	void operator++()
	{
		++determinant;
		++boundaryMass;
		++impedance;
	}

	void move(int n)
	{
		determinant += n;
		boundaryMass += n;
		impedance += n;
	}

	void operator()()
	{
		std::fill(boundaryMass.data, boundaryMass.data + boundaryMass.inc, 0);
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN1M1N<nodes>(determinant[gpindex] * weight[gpindex] / impedance[gpindex], N.data + nodes * gpindex, boundaryMass.data);
		}
	}
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_H_ */
