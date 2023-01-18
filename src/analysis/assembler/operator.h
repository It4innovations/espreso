
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATOR_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATOR_H_

#include "esinfo/eslog.h"
#include "math/simd/simd.h"

namespace espreso {

class Assembler;

struct ActionOperator {
	int isconst, update;

	ActionOperator(): isconst(1), update(1) {}
	virtual ~ActionOperator() {}

	virtual void operator++() {};
	virtual void operator()() {};

	virtual void move(int n)  {};
};

struct InsertOperator {
	virtual ~InsertOperator() {}

	virtual void sisd() =0;
	virtual void simd() =0;
	virtual void peel(size_t size) =0;
	virtual void move(int n)  =0;
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim>
struct HeatTransferOperator {
	virtual ~HeatTransferOperator() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 9];
			alignas(SIMD::size * sizeof(double)) double angle       [SIMD::size * gps * 6];
			alignas(SIMD::size * sizeof(double)) double density     [SIMD::size * gps];
			alignas(SIMD::size * sizeof(double)) double heatCapacity[SIMD::size * gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) double coords[SIMD::size * nodes * ndim];
		alignas(SIMD::size * sizeof(double)) double gpcoords[SIMD::size * gps * ndim];

		alignas(SIMD::size * sizeof(double)) double  w [SIMD::size * gps];
		alignas(SIMD::size * sizeof(double)) double  N [SIMD::size * gps * nodes];
		alignas(SIMD::size * sizeof(double)) double dN [SIMD::size * gps * nodes * edim];

		alignas(SIMD::size * sizeof(double)) double dND[SIMD::size * gps * nodes * edim];
		alignas(SIMD::size * sizeof(double)) double det[SIMD::size * gps];

		alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 9]; // 3 x 3
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

template <size_t nodes, size_t gps, size_t edim>
struct HeatTransferOperator<nodes, gps, 2, edim> {
	virtual ~HeatTransferOperator() {}

	struct Element {
		struct {
			alignas(SIMD::size * sizeof(double)) double thickness[SIMD::size * gps];

			alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 4];
			alignas(SIMD::size * sizeof(double)) double angle       [SIMD::size * gps * 4];
			alignas(SIMD::size * sizeof(double)) double density     [SIMD::size * gps];
			alignas(SIMD::size * sizeof(double)) double heatCapacity[SIMD::size * gps];
		} ecf;

		alignas(SIMD::size * sizeof(double)) double coords[SIMD::size * nodes * 2];
		alignas(SIMD::size * sizeof(double)) double gpcoords[SIMD::size * gps * 2];

		alignas(SIMD::size * sizeof(double)) double  w [SIMD::size * gps];
		alignas(SIMD::size * sizeof(double)) double  N [SIMD::size * gps * nodes];
		alignas(SIMD::size * sizeof(double)) double dN [SIMD::size * gps * nodes * edim];

		alignas(SIMD::size * sizeof(double)) double dND[SIMD::size * gps * nodes * edim];
		alignas(SIMD::size * sizeof(double)) double det[SIMD::size * gps];

		alignas(SIMD::size * sizeof(double)) double conductivity[SIMD::size * gps * 4]; // 2 x 2
	};

	virtual void sisd(Element &element) =0;
	virtual void simd(Element &element) =0;
	virtual void peel(Element &element, size_t size) { simd(element); }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATOR_H_ */
