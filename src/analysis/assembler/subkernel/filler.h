
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_FILLER_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_FILLER_H_

#include "subkernel.h"
#include "math/primitives/matrix_info.h"
#include "math/physics/matrix_base.h"
#include "math/physics/vector_base.h"

namespace espreso {

struct DataFiller: SubKernel {
	const char* name() const { return "MatrixFiller"; }

	size_t dofs, filter, elements;
	double *local, *global;
	const esint *position;
	Matrix_Shape shape;

	DataFiller()
	: dofs(0), filter(0), elements(0),
	  local(nullptr), global(nullptr), position(nullptr),
	  shape(Matrix_Shape::FULL)
	{
		isconst = false;
		action = Assembler::FILL;
	}

	void activate(size_t interval, size_t dofs, size_t elements, double *local, Matrix_Base<double> *A)
	{
		if (A == nullptr) {
			return;
		}
		this->dofs = dofs;
		this->elements = elements;
		this->local = local;
		this->global = A->mapping.elements[interval].data;
		this->position = A->mapping.elements[interval].position;
		this->isactive = 1;
		this->shape = A->shape;
	}

	void activate(size_t interval, size_t dofs, size_t elements, double *local, Vector_Base<double> *A)
	{
		if (A == nullptr) {
			return;
		}
		this->dofs = dofs;
		this->elements = elements;
		this->local = local;
		this->global = A->mapping.elements[interval].data;
		this->position = A->mapping.elements[interval].position;
		this->isactive = 1;
	}

	void activate(size_t region, size_t interval, size_t dofs, size_t elements, double *local, Vector_Base<double> *A)
	{
		if (A == nullptr) {
			return;
		}
		this->dofs = dofs;
		this->filter = A->mapping.boundary[region][interval].filter;
		this->elements = elements;
		this->local = local;
		this->global = A->mapping.boundary[region][interval].data;
		this->position = A->mapping.boundary[region][interval].position;
		this->isactive = 1;
	}
};

template <size_t nodes, class Physics>
struct MatricFillerKernel: DataFiller, Physics {
	MatricFillerKernel(const DataFiller &base): DataFiller(base) {}

	// TODO: direct store
	void simd(typename Physics::Element &element)
	{
		size_t size = nodes * dofs;
		size_t count = std::min((size_t)SIMD::size, elements);
		switch (shape) {
		case Matrix_Shape::FULL:
			for (size_t s = 0, i = 0; s < count; ++s) {
				for (size_t r = 0; r < size; ++r) {
					for (size_t c = 0; c < size; ++c, ++i) {
						global[position[i]] += *(local + (r * size + c) * SIMD::size + s);
					}
				}
			}
			position += count * size * size;
			break;
		case Matrix_Shape::UPPER:
			for (size_t s = 0, i = 0; s < count; ++s) {
				for (size_t r = 0; r < size; ++r) {
					for (size_t c = r; c < size; ++c, ++i) {
						global[position[i]] += *(local + (r * size + c) * SIMD::size + s);
					}
				}
			}
			position += count * ((size * (size - 1)) / 2 + size);
			break;
		case Matrix_Shape::LOWER:
			for (size_t s = 0, i = 0; s < count; ++s) {
				for (size_t r = 0; r < size; ++r) {
					for (size_t c = 0; c <= r; ++c, ++i) {
						global[position[i]] += *(local + (r * size + c) * SIMD::size + s);
					}
				}
			}
			position += count * ((size * (size - 1)) / 2 + size);
			break;
		}
		std::fill(local, local + SIMD::size * size * size, 0.);
		local += count * size * size;
		elements -= count;
	}
};

template <size_t nodes, class Physics>
struct VectorFillerKernel: DataFiller, Physics {
	VectorFillerKernel(const DataFiller &base): DataFiller(base) {}

	void simd(typename Physics::Element &element)
	{
		size_t size = nodes * dofs;
		size_t count = std::min((size_t)SIMD::size, elements);
		for (size_t s = 0, i = 0; s < count; ++s) {
			for (size_t r = 0, j = 0; r < size; ++r, ++i, ++j) {
				global[position[i]] += *(local + j * SIMD::size + s);
			}
		}
		std::fill(local, local + SIMD::size * size, 0.);
		local += count * size;
		position += count * size;
		elements -= count;
	}
};

template <size_t nodes, class Physics>
struct VectorSetterKernel: DataFiller, Physics {
	VectorSetterKernel(const DataFiller &base, const std::function<double(typename Physics::Element &element, size_t &n, size_t &d, size_t &s)> setter): DataFiller(base), setter(setter) {}
	const std::function<double(typename Physics::Element &element, size_t &n, size_t &d, size_t &s)> setter;

	void simd(typename Physics::Element &element)
	{
		size_t size = nodes * dofs;
		size_t count = std::min((size_t)SIMD::size, elements);
		for (size_t s = 0; s < count; ++s) {
			for (size_t d = 0; d < dofs; ++d) {
				if (filter & (1 << d)) {
					for (size_t n = 0; n < nodes; ++n) {
						global[*position++] = setter(element, n, d, s);
					}
				}
			}
		}
		local += count * size;
		elements -= count;
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_FILLER_H_ */
