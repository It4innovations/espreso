
#ifndef SRC_CONFIG_ECF_MATERIAL_TENSOR_H_
#define SRC_CONFIG_ECF_MATERIAL_TENSOR_H_

#include "config/holders/expression.h"
#include "config/description.h"

namespace espreso {

struct TensorConfiguration {

	size_t size;
	std::vector<ECFExpression> values;

	TensorConfiguration(size_t size);

	ECFExpression& get(size_t row, size_t column);
	const ECFExpression& get(size_t row, size_t column) const;

private:
	size_t _get(size_t row, size_t column) const;
};

}

#endif /* SRC_CONFIG_ECF_MATERIAL_TENSOR_H_ */
