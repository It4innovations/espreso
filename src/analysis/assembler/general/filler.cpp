
#include "filler.h"
#include "analysis/math/matrix_base.h"
#include "analysis/math/vector_base.h"

namespace espreso {

void DataFiller::activate(size_t interval, size_t dofs, size_t elements, Matrix_Base<double> *A)
{
    if (A == nullptr) {
        return;
    }
    this->dofs = dofs;
    this->elements = elements;
    this->out = A->mapping.elements[interval].data;
    this->position = A->mapping.elements[interval].position;;
    this->size = A->size();
    this->isactive = 1;
    this->shape = A->shape;
}

void DataFiller::activate(size_t interval, size_t dofs, size_t elements, Vector_Base<double> *A)
{
    if (A == nullptr) {
        return;
    }
    this->dofs = dofs;
    this->elements = elements;
    this->out = A->mapping.elements[interval].data;
    this->position = A->mapping.elements[interval].position;
    this->size = A->size();
    this->isactive = 1;
}

void DataFiller::activate(size_t region, size_t interval, size_t dofs, size_t elements, Vector_Base<double> *A)
{
    if (A == nullptr) {
        return;
    }
    this->dofs = dofs;
    this->filter = A->mapping.boundary[region][interval].filter;
    this->elements = elements;
    this->out = A->mapping.boundary[region][interval].data;
    this->position = A->mapping.boundary[region][interval].position;
    this->size = A->size();
    this->isactive = 1;
}

}
