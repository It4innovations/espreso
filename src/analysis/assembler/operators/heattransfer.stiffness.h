
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

	void move(size_t n)
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

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness2<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC, Physics>: HeatTransferStiffnessBase, Physics {
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

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness2<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferStiffnessBase, Physics {
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

template <size_t nodes, size_t gps, size_t edim, class Physics>
struct HeatTransferStiffness2<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, Physics>: HeatTransferStiffnessBase, Physics {
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
