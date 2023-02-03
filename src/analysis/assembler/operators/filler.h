
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_FILLER_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "math/physics/matrix_base.h"

namespace espreso {

struct MatrixFiller: ActionOperator {
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

	void sisd(typename Physics::Element &element)
	{
		for (size_t r = 0, i = 0; r < nodes * dofs; ++r) {
			for (size_t c = r; c < nodes * dofs; ++c, ++i) {
				global[position[i]] += *(local.data + r * nodes * dofs + c);
			}
		}
		std::fill(local.data, local.data + local.inc, 0.);
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0, i = 0; s < SIMD::size; ++s) {
			for (size_t r = 0; r < nodes * dofs; ++r) {
				for (size_t c = r; c < nodes * dofs; ++c, ++i) {
					global[position[i]] += *(local.data + (r * nodes * dofs + c) * SIMD::size + s);
				}
			}
		}
		std::fill(local.data, local.data + SIMD::size * local.inc, 0.);
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0, i = 0; s < size; ++s) {
			for (size_t r = 0; r < nodes * dofs; ++r) {
				for (size_t c = r; c < nodes * dofs; ++c, ++i) {
					global[position[i]] += *(local.data + (r * nodes * dofs + c) * SIMD::size + s);
				}
			}
		}
		std::fill(local.data, local.data + size * local.inc, 0.);
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

	void sisd(typename Physics::Element &element)
	{
		for (size_t r = 0, i = 0; r < nodes * dofs; ++r) {
			for (size_t c = 0; c < nodes * dofs; ++c, ++i) {
				global[position[i]] += *(local.data + r * nodes * dofs + c);
			}
		}
		std::fill(local.data, local.data + local.inc, 0.);
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0, i = 0; s < SIMD::size; ++s) {
			for (size_t r = 0; r < nodes * dofs; ++r) {
				for (size_t c = 0; c < nodes * dofs; ++c, ++i) {
					global[position[i]] += *(local.data + (r * nodes * dofs + c) * SIMD::size + s);
				}
			}
		}
		std::fill(local.data, local.data + SIMD::size * local.inc, 0.);
		move(SIMD::size);
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
		std::fill(local.data, local.data + size * local.inc, 0.);
		move(size);
	}
};

template<size_t nodes, size_t gps, size_t dimension>
struct VectorFiller: ActionOperator {
	VectorFiller(int interval, const ParameterData &rhs, double *global, const esint *position)
	: rhs(rhs, interval),
	  global(global), position(position) {}

	InputParameterIterator rhs;
	double *global;
	const esint *position;

	void operator++()
	{
		++rhs;
	}

	void move(int n)
	{
		rhs += n;
		position += n * nodes * dimension;
	}

	void operator()()
	{
		for (size_t r = 0; r < nodes * dimension; ++r) {
			global[*position++] += rhs[r];
		}
	}
};

// Setter never rewrite output since output is directly rewriten by expression
template <size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype, class Physics, class Setter>
struct VectorSetter: ActionOperator, Physics {
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
		action = ActionOperator::Action::ASSEMBLE;
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

	void sisd(typename Physics::Element &element)
	{
		for (size_t d = 0, i = 0; d < dofs; ++d) {
			for (size_t n = 0; n < nodes; ++n, ++i) {
				if (filter & (1 << d)) {
					global[position[i]] = setter(element, n, d, 0);
				}
			}
		}
		move(1);
	}

	void simd(typename Physics::Element &element)
	{
		for (size_t s = 0, i = 0; s < SIMD::size; ++s) {
			for (size_t d = 0; d < dofs; ++d) {
				for (size_t n = 0; n < nodes; ++n, ++i) {
					if (filter & (1 << d)) {
						global[position[i]] = setter(element, n, d, s);
					}
				}
			}
		}
		move(SIMD::size);
	}

	void peel(typename Physics::Element &element, size_t size)
	{
		for (size_t s = 0, i = 0; s < size; ++s) {
			for (size_t d = 0; d < dofs; ++d) {
				for (size_t n = 0; n < nodes; ++n, ++i) {
					if (filter & (1 << d)) {
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
