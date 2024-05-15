
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_FILLER_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_FILLER_H_

#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"
#include "subkernel.h"
#include "math/primitives/matrix_info.h"
#include <functional>

namespace espreso {

struct DataFiller: SubKernel {
    const char* name() const { return "MatrixFiller"; }

    size_t dofs, filter, elements;
    double *out;
    const esint *position;
    Matrix_Shape shape;

    DataFiller()
    : dofs(0), filter(0), elements(0),
      out(nullptr), position(nullptr),
      shape(Matrix_Shape::FULL)
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(size_t interval, size_t dofs, size_t elements, Matrix_Base<double> *A)
    {
        if (A == nullptr) {
            return;
        }
        this->dofs = dofs;
        this->elements = elements;
        this->out = A->mapping.elements[interval].data;
        this->position = A->mapping.elements[interval].position;
        this->isactive = 1;
        this->shape = A->shape;
    }

    void activate(size_t interval, size_t dofs, size_t elements, Vector_Base<double> *A)
    {
        if (A == nullptr) {
            return;
        }
        this->dofs = dofs;
        this->elements = elements;
        this->out = A->mapping.elements[interval].data;
        this->position = A->mapping.elements[interval].position;
        this->isactive = 1;
    }

    void activate(size_t region, size_t interval, size_t dofs, size_t elements, Vector_Base<double> *A)
    {
        if (A == nullptr) {
            return;
        }
        this->dofs = dofs;
        this->filter = A->mapping.boundary[region][interval].filter;
        this->elements = elements;
        this->out = A->mapping.boundary[region][interval].data;
        this->position = A->mapping.boundary[region][interval].position;
        this->isactive = 1;
    }
};

template <size_t nodes>
struct MatrixFillerKernel: DataFiller {
    MatrixFillerKernel(const DataFiller &base): DataFiller(base) {}

    void simd(SIMD matrix[])
    {
        size_t size = nodes * dofs;
        size_t count = std::min((size_t)SIMD::size, elements);
        switch (shape) {
        case Matrix_Shape::FULL:
            for (size_t s = 0, i = 0; s < count; ++s) {
                for (size_t r = 0; r < size; ++r) {
                    for (size_t c = 0; c < size; ++c, ++i) {
                        out[position[i]] += matrix[(r * size + c)][s];
                    }
                }
            }
            position += count * size * size;
            break;
        case Matrix_Shape::UPPER:
            for (size_t s = 0, i = 0; s < count; ++s) {
                for (size_t r = 0; r < size; ++r) {
                    for (size_t c = r; c < size; ++c, ++i) {
                        out[position[i]] += matrix[(r * size + c)][s];
                    }
                }
            }
            position += count * ((size * (size - 1)) / 2 + size);
            break;
        case Matrix_Shape::LOWER:
            for (size_t s = 0, i = 0; s < count; ++s) {
                for (size_t r = 0; r < size; ++r) {
                    for (size_t c = 0; c <= r; ++c, ++i) {
                        out[position[i]] += matrix[(r * size + c)][s];
                    }
                }
            }
            position += count * ((size * (size - 1)) / 2 + size);
            break;
        }
        elements -= count;
        for (size_t i = 0; i < size * size; ++i) {
            matrix[i] = zeros();
        }
    }
};

template <size_t nodes>
struct RHSFillerKernel: DataFiller {
    RHSFillerKernel(const DataFiller &base): DataFiller(base) {}

    void simd(SIMD vector[])
    {
        size_t size = nodes * dofs;
        size_t count = std::min((size_t)SIMD::size, elements);
        for (size_t s = 0, i = 0; s < count; ++s) {
            for (size_t r = 0; r < size; ++r, ++i) {
                out[position[i]] += vector[r][s];
            }
        }
        position += count * size;
        elements -= count;
        for (size_t i = 0; i < size; ++i) {
            vector[i] = zeros();
        }
    }
};

template <size_t nodes, class Element>
struct VectorSetterKernel: DataFiller, Element {
    VectorSetterKernel(const DataFiller &base, const std::function<double(Element &element, size_t &n, size_t &d, size_t &s)> setter): DataFiller(base), setter(setter) {}
    const std::function<double(Element &element, size_t &n, size_t &d, size_t &s)> setter;

    void simd(Element &element)
    {
        size_t count = std::min((size_t)SIMD::size, elements);
        for (size_t s = 0; s < count; ++s) {
            for (size_t d = 0; d < dofs; ++d) {
                if (filter & (1 << d)) {
                    for (size_t n = 0; n < nodes; ++n) {
                        out[*position++] = setter(element, n, d, s);
                    }
                }
            }
        }
        elements -= count;
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_FILLER_H_ */
