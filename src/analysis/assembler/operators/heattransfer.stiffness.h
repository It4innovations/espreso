
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"
#include "analysis/assembler/module/heattransfer.element.h"

namespace espreso {

struct HeatTransferStiffnessBase: public ActionOperator {
	OutputParameterIterator stiffness;

	HeatTransferStiffnessBase(size_t interval, ParameterData &stiffness)
	: stiffness(stiffness, interval)
	{
		isconst = false;
	}

	void move(int n)
	{
		stiffness += n;
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics> struct HeatTransferStiffness2;

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness2<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			double scale = element.ecf.thickness[gpindex][0] * element.det[gpindex][0] * element.w[gpindex][0] * element.conductivity[gpindex][0][0];
			double scale = element.det[gpindex] * element.w[gpindex];
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t m = 0; m < nodes; ++m) {
					stiffness[n * nodes + m] += scale * (element.dND[gpindex][n][0] * element.dND[gpindex][m][0] + element.dND[gpindex][n][1] * element.dND[gpindex][m][1]);
				}
			}
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			SIMD scale = load(element.ecf.thickness[gpindex]) * load(element.det[gpindex]) * load(element.w[gpindex]) * load(element.conductivity[gpindex][0]);
			SIMD scale = element.det[gpindex] * element.w[gpindex];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gpindex][n][0];
				SIMD ny = element.dND[gpindex][n][1];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gpindex][m][0];
					SIMD my = element.dND[gpindex][m][1];
					SIMD res = load(stiffness.data + (n * nodes + m) * SIMD::size);
					res = res + scale * (nx * mx + ny * my);
					store(stiffness.data + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness2<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			double scale = element.det[gpindex][0] * element.w[gpindex][0] * element.conductivity[gpindex][0][0];
			double scale = element.det[gpindex] * element.w[gpindex];
			for (size_t n = 0; n < nodes; ++n) {
				for (size_t m = 0; m < nodes; ++m) {
					stiffness[n * nodes + m] += scale * (
							element.dND[gpindex][n][0] * element.dND[gpindex][m][0] +
							element.dND[gpindex][n][1] * element.dND[gpindex][m][1] +
							element.dND[gpindex][n][2] * element.dND[gpindex][m][2]);
				}
			}
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		double * __restrict__ out = stiffness.data;
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			SIMD scale = load(element.det[gpindex]) * load(element.w[gpindex]) * load(element.conductivity[gpindex]);
			SIMD scale = element.det[gpindex] * element.w[gpindex];
			for (size_t n = 0; n < nodes; ++n) {
				SIMD nx = element.dND[gpindex][n][0];
				SIMD ny = element.dND[gpindex][n][1];
				SIMD nz = element.dND[gpindex][n][2];
				for (size_t m = 0; m < nodes; ++m) {
					SIMD mx = element.dND[gpindex][m][0];
					SIMD my = element.dND[gpindex][m][1];
					SIMD mz = element.dND[gpindex][m][2];
					SIMD res = load(out + (n * nodes + m) * SIMD::size);

					res = res + scale * (nx * mx + ny * my + nz * mz);
					store(out + (n * nodes + m) * SIMD::size, res);
				}
			}
		}
		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness2<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDMN2M22M2N<nodes>(element.ecf.thickness[gpindex] * element.det[gpindex] * element.w[gpindex], element.conductivity + 4 * gpindex, element.dND + 2 * nodes * gpindex, stiffness.data);
//		}
//		move(1);
	}

	void simd(typename Physics::Element &element)
	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			SIMD scale =  load(element.ecf.thickness + gpindex * SIMD::size)
//						* load(element.det + gpindex * SIMD::size)
//						* load(element.w + gpindex * SIMD::size);
//
//			ADDMN2M22M2NSimd<nodes>(scale, element.conductivity + 4 * gpindex * SIMD::size, element.dND + 2 * nodes * gpindex * SIMD::size, stiffness.data);
//		}
//		move(SIMD::size);
	}
};

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness2<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			ADDMN3M33M3N<nodes>(element.det[gpindex] * element.w[gpindex], element.conductivity + 9 * gpindex, element.dND + 3 * nodes * gpindex, stiffness.data);
//		}
//		move(1);
	}

	void simd(typename Physics::Element &element)
	{
//		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			SIMD scale =  load(element.det + gpindex * SIMD::size)
//						* load(element.w + gpindex * SIMD::size);
//
//			ADDMN3M33M3NSimd<nodes>(scale, element.conductivity + 9 * gpindex * SIMD::size, element.dND + 3 * nodes * gpindex * SIMD::size, stiffness.data);
//		}
//		move(SIMD::size);
	}
};

template<size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct HeatTransferStiffnessIsotropic;

template<size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffnessIsotropic<nodes, gps, 2, edim, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN2M2N<nodes>(element.ecf.thickness[gpindex] * element.det[gpindex] * element.w[gpindex] * element.conductivity[gpindex], element.dND + 2 * nodes * gpindex, stiffness.data);
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			SIMD scale =  load(element.ecf.thickness + gpindex * SIMD::size)
//						* load(element.det + gpindex * SIMD::size)
//						* load(element.w + gpindex * SIMD::size)
//						* load(element.conductivity + gpindex * SIMD::size);
			SIMD scale =  load(element.det + gpindex * SIMD::size)
						* load(element.w + gpindex * SIMD::size);
			ADDMN2M2NSimd<nodes>(scale, element.dND + 2 * nodes * gpindex * SIMD::size, stiffness.data);
		}
		move(SIMD::size);
	}
};

template<size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffnessIsotropic<nodes, gps, 3, edim, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN3M3N<nodes>(element.det[gpindex] * element.w[gpindex] * element.conductivity[gpindex], element.dND + 3 * nodes * gpindex, stiffness.data);
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
//			SIMD scale =  load(element.det + gpindex * SIMD::size)
//						* load(element.w + gpindex * SIMD::size)
//						* load(element.conductivity + gpindex * SIMD::size);
			SIMD scale =  load(element.det + gpindex * SIMD::size)
						* load(element.w + gpindex * SIMD::size);
			ADDMN3M3NSimd<nodes>(scale, element.dND + 3 * nodes * gpindex * SIMD::size, stiffness.data);
		}
		move(SIMD::size);
	}
};

template<size_t nodes, size_t gps, size_t ndim, size_t edim, class Physics> struct HeatTransferStiffness;

template<size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness<nodes, gps, 2, edim, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN2M22M2N<nodes>(element.ecf.thickness[gpindex] * element.det[gpindex] * element.w[gpindex], element.conductivity + 4 * gpindex, element.dND + 2 * nodes * gpindex, stiffness.data);
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			SIMD scale =  load(element.ecf.thickness + gpindex * SIMD::size)
						* load(element.det + gpindex * SIMD::size)
						* load(element.w + gpindex * SIMD::size);

			ADDMN2M22M2NSimd<nodes>(scale, element.conductivity + 4 * gpindex * SIMD::size, element.dND + 2 * nodes * gpindex * SIMD::size, stiffness.data);
		}
		move(SIMD::size);
	}
};

template<size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness<nodes, gps, 3, edim, Physics>: HeatTransferStiffnessBase, Physics {
	using HeatTransferStiffnessBase::HeatTransferStiffnessBase;

	void sisd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			ADDMN3M33M3N<nodes>(element.det[gpindex] * element.w[gpindex], element.conductivity + 9 * gpindex, element.dND + 3 * nodes * gpindex, stiffness.data);
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			SIMD scale =  load(element.det + gpindex * SIMD::size)
						* load(element.w + gpindex * SIMD::size);

			ADDMN3M33M3NSimd<nodes>(scale, element.conductivity + 9 * gpindex * SIMD::size, element.dND + 3 * nodes * gpindex * SIMD::size, stiffness.data);
		}
		move(SIMD::size);
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_HEATTRANSFER_STIFFNESS_H_ */
