
#include "tensor.h"

using namespace espreso;

TensorConfiguration::TensorConfiguration(size_t size)
: size(size), values(size * size)
{

}

size_t TensorConfiguration::_get(size_t row, size_t column) const
{
	return row * size + column;
}

const ECFExpression& TensorConfiguration::get(size_t row, size_t column) const
{
	return values[_get(row, column)];
}

ECFExpression& TensorConfiguration::get(size_t row, size_t column)
{
	return values[_get(row, column)];
}



