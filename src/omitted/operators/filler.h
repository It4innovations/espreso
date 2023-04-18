
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include <analysis/assembler/subkernel/operator.h>
#include "analysis/assembler/parameter.h"
#include "math/simd/simd.h"
#include "math/physics/matrix_base.h"

namespace espreso {

struct MatrixFiller: ActionOperator {
	const char* name() const { return "MatrixFiller"; }

	MatrixFiller(size_t interval, size_t dofs, ParameterData &local, Matrix_Base<double> *A)
	: dofs(dofs),
	  local(local, interval),
	  global(A->mapping.elements[interval].data), position(A->mapping.elements[interval].position)
	{
		isconst = false;
		action = ActionOperator::Action::FILL;
	}

	size_t dofs;
	OutputParameterIterator local;
	double *global;
	const esint *position;
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct SymmetricMatricFiller: MatrixFiller, Physics {
	using MatrixFiller::MatrixFiller;

	void move(int n)
	{
		local += n;
		position += n * ((nodes * dofs * (nodes * dofs - 1)) / 2 + nodes * dofs);
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0, i = 0; s < size; ++s) {
			for (size_t r = 0, j = 0; r < nodes * dofs; ++r) {
				for (size_t c = r; c < nodes * dofs; ++c, ++i, ++j) {
					global[position[i]] += *(local.data + j * SIMD::size + s);
				}
			}
		}
		std::fill(local.data, local.data + SIMD::size * local.inc, 0.);
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct SymmetricToFullMatricFiller: MatrixFiller, Physics {
	using MatrixFiller::MatrixFiller;

	void move(int n)
	{
		local += n;
		position += n * nodes * dofs * nodes * dofs;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0, i = 0; s < size; ++s) {
			for (size_t r = 0; r < nodes * dofs; ++r) {
				for (size_t c = 0; c < r; ++c, ++i) {
					global[position[i]] += *(local.data + (c * nodes * dofs + r - ((c + 1) * c / 2)) * SIMD::size + s);
				}
				for (size_t c = r; c < nodes * dofs; ++c, ++i) {
					global[position[i]] += *(local.data + (r * nodes * dofs + c - ((r + 1) * r / 2)) * SIMD::size + s);
				}
			}
		}
		std::fill(local.data, local.data + SIMD::size * local.inc, 0.);
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct GeneralMatricFiller: MatrixFiller, Physics {
	using MatrixFiller::MatrixFiller;

	void move(int n)
	{
		local += n;
		position += n * nodes * dofs * nodes * dofs;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0, i = 0; s < size; ++s) {
			for (size_t r = 0; r < nodes * dofs; ++r) {
				for (size_t c = 0; c < nodes * dofs; ++c, ++i) {
					global[position[i]] += *(local.data + (r * nodes * dofs + c) * SIMD::size + s);
				}
			}
		}
		std::fill(local.data, local.data + SIMD::size * local.inc, 0.);
		move(size);
	}
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics>
struct VectorFiller: ActionOperator, Physics {
	const char* name() const { return "VectorFiller"; }

	VectorFiller(size_t interval, size_t dofs, ParameterData &local, Vector_Base<double> *v)
	: dofs(dofs),
	  local(local, interval),
	  global(v->mapping.elements[interval].data),
	  position(v->mapping.elements[interval].position)
	{
		isconst = false;
		action = ActionOperator::Action::FILL;
	}

	VectorFiller(size_t region, size_t interval, size_t dofs, ParameterData &local, Vector_Base<double> *v)
	: dofs(dofs),
	  local(local, interval),
	  global(v->mapping.boundary[region][interval].data),
	  position(v->mapping.boundary[region][interval].position)
	{
		isconst = false;
		action = ActionOperator::Action::FILL;
	}

	size_t dofs;
	OutputParameterIterator local;
	double *global;
	const esint *position;

	void move(int n)
	{
		local += n;
		position += n * nodes * dofs;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0, i = 0; s < size; ++s) {
			for (size_t d = 0; d < dofs; ++d) {
				for (size_t n = 0; n < nodes; ++n, ++i) {
					global[position[i]] += *(local.data + (d * nodes + n) * SIMD::size + s);
				}
			}
		}
		std::fill(local.data, local.data + SIMD::size * local.inc, 0.);
		move(size);
	}
};

// Setter never rewrite output since output is directly rewriten by expression
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct VectorSetter: ActionOperator, Physics {
	const char* name() const { return "VectorSetter"; }

	VectorSetter(size_t region, size_t interval, size_t dofs, Vector_Base<double> *v, const Setter &setter)
	: dofs(dofs), dim(0),
	  global(v->mapping.boundary[region][interval].data),
	  position(v->mapping.boundary[region][interval].position),
	  filter(v->mapping.boundary[region][interval].filter),
	  setter(setter)
	{
		for (size_t d = 0; d < dofs; ++d) {
			if (filter & (1 << d)) {
				++dim;
			}
		}
		isconst = false;
		action = ActionOperator::Action::ASSEMBLE | Action::REASSEMBLE;
	}

	size_t dofs, dim;
	double *global;
	const esint *position;
	const int filter;
	const Setter setter;

	void move(int n)
	{
		position += n * nodes * dim;
	}

	void simd(typename Physics::Element &element)
	{
		peel(element, SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0, i = 0; s < size; ++s) {
			for (size_t d = 0; d < dofs; ++d) {
				if (filter & (1 << d)) {
					for (size_t n = 0; n < nodes; ++n, ++i) {
						global[position[i]] = setter(element, n, d, s);
					}
				}
			}
		}
		move(size);
	}
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_ */
