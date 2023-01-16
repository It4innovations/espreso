
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_FORCES_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_FORCES_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

template<size_t nodes, size_t gps>
struct HeatRHS: public ActionOperator {
	HeatRHS(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &heatSource, ParameterData &rhs)
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
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[n] += J[gpindex] * weight[gpindex] * heatSource[gpindex] * N[gpindex * nodes + n];
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

template<size_t nodes, size_t gps>
struct HeatQ: public ActionOperator {
	HeatQ(int interval, double area, const ParameterData &heatFlow, const ParameterData &heatFlux, const ParameterData &htc, const ParameterData &extTemp, ParameterData &q)
	: area(area), heatFlow(heatFlow, interval),
	  heatFlux(heatFlux, interval),
	  htc(htc, interval), extTemp(extTemp, interval),
	  q(q, interval)
	{

	}

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			q[gpindex] += heatFlow[gpindex] / area;
			q[gpindex] += heatFlux[gpindex];
			q[gpindex] += htc[gpindex] * extTemp[gpindex];
		}
	}

	void operator++()
	{
		++heatFlow; ++heatFlux;
		++htc; ++extTemp;
		++q;
	}

	void move(int n)
	{
		heatFlow += n; heatFlux += n;
		htc += n; extTemp += n;
		q += n;
	}

	double area;
	InputParameterIterator heatFlow, heatFlux, htc, extTemp;
	OutputParameterIterator q;
};

template<size_t nodes, size_t gps>
struct HeatRHS2D: public ActionOperator {
	HeatRHS2D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &thickness, const ParameterData &q, ParameterData &rhs)
	: N(N, interval),
	  weight(weight, interval),
	  J(J, interval),
//	  thickness(thickness, interval), // TODO
	  q(q, interval),
	  rhs(rhs, interval)
	{ }

	InputParameterIterator N, weight, J;
//	InputParameterIterator thickness;
	InputParameterIterator q;
	OutputParameterIterator rhs;

	void operator()()
	{
		for (size_t n = 0; n < nodes; ++n) {
			for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
				rhs[n] += J[gpindex] * weight[gpindex] * q[gpindex] * N[gpindex * nodes + n];
//				rhs[n] += thickness[gpindex] * J[gpindex] * weight[gpindex] * q[gpindex] * N[gpindex * nodes + n];
			}
		}
	}

	void operator++()
	{
		++weight; ++J;
//		++thickness;
		++q;
		++rhs;
	}

	void move(int n)
	{
		weight += n; J += n;
//		thickness += n;
		q += n;
		rhs += n;
	}
};

template<size_t nodes, size_t gps>
struct HeatRHS3D: public ActionOperator {
	HeatRHS3D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &q, ParameterData &rhs)
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

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_FORCES_H_ */
