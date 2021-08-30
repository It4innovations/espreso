
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
	: dND(dND, interval),
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

	void move(int n)
	{
		dND += n; determinant += n;
		stiffness += n;
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
	: N(N, interval),
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

	void move(int n)
	{
		N += n; determinant += n;
		mass += n;
	}

	void operator()()
	{
		// double kappa = omega / speedOfSound
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN1M1N<nodes>(determinant[gpindex] * weight[gpindex] /*  * kappa[gpindex] * kappa[gpindex] */, N.data + nodes * gpindex, mass.data);
		}
	}
};

template <size_t nodes, size_t gps>
struct AcousticQ: public ActionOperator {
	AcousticQ(int interval, double area, const ParameterData &g,  ParameterData &impedance, ParameterData &q)
	: ActionOperator(interval, g.isconst[interval], g.update[interval]),
	  area(area),
	  g(g, interval),
	  impedance(impedance, interval),
	  q(q, interval)
	{
//		if (update) {
//			std::fill((q.data->begin() + interval)->data(), (q.data->begin() + interval + 1)->data(), 0);
//		}
	}

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			q.data[gpindex] += g.data[gpindex] / area;
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
//		if (update) {
//			std::fill((rhs.data->begin() + interval)->data(), (rhs.data->begin() + interval + 1)->data(), 0);
//		}
	}

	InputParameterIterator N, weight, J;
	InputParameterIterator q;
	OutputParameterIterator rhs;

	void operator()()
	{
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
	: ActionOperator(interval, boundaryMass.isconst[interval], boundaryMass.update[interval]),
	  N(N, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  impedance(impedance, interval),
	  boundaryMass(boundaryMass, interval)
	{
//		if (update) {
//			std::fill((boundaryMass.data->begin() + interval)->data(), (boundaryMass.data->begin() + interval + 1)->data(), 0);
//		}
	}

	void operator++()
	{
		++impedance, ++N; ++determinant;
		++boundaryMass;
	}

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN1M1N<nodes>(determinant[gpindex] * weight[gpindex], N.data + nodes * gpindex, boundaryMass.data);
		}
	}

	void reset()
	{

	}
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_H_ */
