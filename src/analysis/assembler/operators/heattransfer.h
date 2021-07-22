
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_H_

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

template<size_t nodes, size_t gps>
struct HeatQ: public ActionOperator {
	HeatQ(int interval, double area, const ParameterData &heatFlow, const ParameterData &heatFlux, const ParameterData &htc, const ParameterData &extTemp, ParameterData &q)
	: ActionOperator(interval, q.isconst[interval], q.update[interval]),
	  area(area), heatFlow(heatFlow, interval),
	  heatFlux(heatFlux, interval),
	  htc(htc, interval), extTemp(extTemp, interval),
	  q(q, interval)
	{
//		if (update) {
//			std::fill((q.data->begin() + interval)->data(), (q.data->begin() + interval + 1)->data(), 0);
//		}
	}

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			q.data[gpindex] += heatFlow.data[gpindex] / area;
			q.data[gpindex] += heatFlux.data[gpindex];
			q.data[gpindex] += htc.data[gpindex] * extTemp.data[gpindex];
			printf("q[%lu] = %f\n", gpindex, q.data[gpindex]);
		}
	}

	void operator++()
	{
		++heatFlow; ++heatFlux;
		++htc; ++extTemp;
		++q;
	}

	void reset()
	{

	}

	double area;
	InputParameterIterator heatFlow, heatFlux, htc, extTemp;
	OutputParameterIterator q;
};

template<size_t nodes, size_t gps>
struct HeatRHS2D: public ActionOperator {
	HeatRHS2D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &thickness, const ParameterData &q, ParameterData &rhs)
	: ActionOperator(interval, rhs.isconst[interval], rhs.update[interval]),
	  N(N, interval),
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
				printf("J=%f w=%f q=%f N=%f\n", J.data[gpindex], weight.data[gpindex], q.data[gpindex], N.data[gpindex * nodes + n]);
				rhs.data[n] += J.data[gpindex] * weight.data[gpindex] * q.data[gpindex] * N.data[gpindex * nodes + n];
//				rhs.data[n] += thickness.data[gpindex] * J.data[gpindex] * weight.data[gpindex] * q.data[gpindex] * N.data[gpindex * nodes + n];
				printf("rhs[%lu] = %f\n", n, rhs.data[n]);
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

	void reset()
	{

	}
};

template<size_t nodes, size_t gps>
struct HeatRHS3D: public ActionOperator {
	HeatRHS3D(int interval, const ParameterData &N, const ParameterData &weight, const ParameterData &J, const ParameterData &q, ParameterData &rhs)
	: ActionOperator(interval, rhs.isconst[interval], rhs.update[interval]),
	  N(N, interval),
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

	void reset()
	{

	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_H_ */
