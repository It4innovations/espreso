
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

namespace espreso {

struct HeatStiffness: public ActionOperator {
	HeatStiffness(
			int interval,
			const ParameterData &dND,
			const ParameterData &weight,
			const ParameterData &determinant,
			const ParameterData &conductivity,
			const ParameterData &xi,
			const ParameterData &thickness,
			ParameterData &stiffness)
	: dND(dND, interval),
	  weight(weight, interval, 0),
	  determinant(determinant, interval),
	  conductivity(conductivity, interval),
	  xi(xi, interval),
	  thickness(thickness, interval),
	  stiffness(stiffness, interval)
	{

	}

	InputParameterIterator dND, weight, determinant, conductivity, xi, thickness;
	OutputParameterIterator stiffness;

	void operator++()
	{
		++dND; ++determinant; ++conductivity; ++xi; ++thickness;
		++stiffness;
	}

	void move(int n)
	{
		dND += n; determinant += n; conductivity += n; xi += n; thickness += n;
		stiffness += n;
	}

	HeatStiffness& operator+=(const size_t rhs)
	{
		dND += rhs; determinant += rhs; conductivity += rhs; xi += rhs; thickness += rhs;
		stiffness += rhs;
		return *this;
	}
};

template<size_t nodes, size_t gps>
struct Stiffness2DHeatIsotropic: public HeatStiffness {
	using HeatStiffness::HeatStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN2M2N<nodes>(thickness[gpindex] * xi[gpindex] * determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}

	void simd()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			SIMD scale =  load(&thickness[gpindex * SIMD::size])
						* load(&determinant[gpindex * SIMD::size])
						* load(&weight[gpindex * SIMD::size])
						* load(&conductivity[gpindex * SIMD::size]);
			ADDMN2M2NSimd<nodes>(scale, dND.data + 2 * nodes * gpindex * SIMD::size, stiffness.data);
		}
		move(SIMD::size);
	}
};

template<size_t nodes, size_t gps>
struct Stiffness2DHeat: public HeatStiffness {
	using HeatStiffness::HeatStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN2M22M2N<nodes>(thickness[gpindex] * xi[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 4 * gpindex, dND.data + 2 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness3DHeatIsotropic: public HeatStiffness {
	using HeatStiffness::HeatStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN3M3N<nodes>(xi[gpindex] * determinant[gpindex] * weight[gpindex] * conductivity[gpindex], dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

template<size_t nodes, size_t gps>
struct Stiffness3DHeat: public HeatStiffness {
	using HeatStiffness::HeatStiffness;

	void operator()()
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN3M33M3N<nodes>(xi[gpindex] * determinant[gpindex] * weight[gpindex], conductivity.data + 9 * gpindex, dND.data + 3 * nodes * gpindex, stiffness.data);
		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_ */
