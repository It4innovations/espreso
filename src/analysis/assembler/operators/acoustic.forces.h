
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_FORCES_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_FORCES_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"

#include <iostream>
#include <complex>

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

		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDM12M2N<nodes>(determinant[gpindex] * weight[gpindex] / density[gpindex], q.data + 2 * gpindex, dND.data + gpindex * nodes * 2, dipole.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct AcousticDipole3D: public AcousticDipole {
	using AcousticDipole::AcousticDipole;
	
	void operator()()
	{
		std::fill(dipole.data, dipole.data + dipole.inc, 0);

		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDM13M3N<nodes>(determinant[gpindex] * weight[gpindex] / density[gpindex], q.data + 3 * gpindex, dND.data + gpindex * nodes * 3, dipole.data);
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
			q[gpindex] += g[gpindex];
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
struct AcousticAcceleration2D: public ActionOperator {

	InputParameterIterator N, weight, J;
	InputParameterIterator normals, acceleration_vector;
	OutputParameterIterator rhs;

	AcousticAcceleration2D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &normals, const ParameterData &acceleration_vector, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  normals(normals, interval),
	  acceleration_vector(acceleration_vector, interval),
	  rhs(rhs, interval)
	{ }

	void operator++()
	{
		++weight; ++J;
		++normals;
		++acceleration_vector;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n;
		normals += n;
		acceleration_vector += n;
		rhs += n;
	}

	void operator()()
	{
		std::fill(rhs.data, rhs.data + rhs.inc, 0);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				double proj = normals[2*gpindex + 0] * acceleration_vector[2*gpindex + 0] + normals[2*gpindex + 1] * acceleration_vector[2*gpindex + 1];
				rhs[n] += J[gpindex] * weight[gpindex] * (-proj) * N[gpindex * nodes + n];
			}
		}
	}
};

template <size_t nodes, size_t gps>
struct AcousticAcceleration3D: public ActionOperator {

	InputParameterIterator N, weight, J;
	InputParameterIterator normals, acceleration_vector;
	OutputParameterIterator rhs;

	AcousticAcceleration3D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &normals, const ParameterData &acceleration_vector, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
	  normals(normals, interval),
	  acceleration_vector(acceleration_vector, interval),
	  rhs(rhs, interval)
	{ }

	void operator++()
	{
		++weight; ++J;
		++normals;
		++acceleration_vector;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n;
		normals += n;
		acceleration_vector += n;
		rhs += n;
	}

	void operator()()
	{
		std::fill(rhs.data, rhs.data + rhs.inc, 0);
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				double proj = normals[3*gpindex + 0] * acceleration_vector[3*gpindex + 0] + normals[3*gpindex + 1] * acceleration_vector[3*gpindex + 1]  + normals[3*gpindex + 2] * acceleration_vector[3*gpindex + 2];
				rhs[n] += J[gpindex] * weight[gpindex] * (-proj) * N[gpindex * nodes + n];
			}
		}
	}
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
				rhs[n] += J[gpindex] * weight[gpindex] * q[gpindex] * N[gpindex * nodes + n];
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
				rhs[n] += J[gpindex] * weight[gpindex] * q[gpindex] * N[gpindex * nodes + n];
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
//			ADDMN1M1N<nodes>(determinant[gpindex] * weight[gpindex] / impedance[gpindex], N.data + nodes * gpindex, boundaryMass.data);
		}
	}
};

template <size_t nodes>
struct AcousticsPointSource_Flow: public ActionOperator {
	InputParameterIterator flowRate, phase;
	OutputParameterIterator sourceAmplitude;

	AcousticsPointSource_Flow(
		int interval,
		const ParameterData &flowRate,
		const ParameterData &phase,
		ParameterData &sourceAmplitude)
	: flowRate(flowRate, interval),
	  phase(phase, interval),
	  sourceAmplitude(sourceAmplitude, interval)
	{

	}

	void operator++()
	{
		++flowRate;
		++phase;
		++sourceAmplitude;
	}

	void move(int n)
	{
		flowRate += n;
		phase += n;
		sourceAmplitude += n;
	}

	void operator()()
	{
		const double omega = 100;
		const double density = 1.25;
		const std::complex<double> i(0, 1);

		for (size_t n = 0; n < nodes; ++n) {
			std::complex<double> amplitude = std::exp(i * phase[n]) * i * omega * density *  flowRate[n] / (4.0 * M_PI);
			this->rhs[n] = amplitude.real();
			// ?? pridat imag slozku
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_ACOUSTIC_FORCES_H_ */
