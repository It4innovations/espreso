
#include "linear.h"

namespace assembler {

template <>
void Linear<API>::KMf(size_t part, bool dynamics)
{
	_K[part] = *(this->_input.K);
}

template <>
void Linear<API>::RHS()
{
	_f[0] = std::vector<double>(this->_input.rhs, this->_input.rhs + this->_input.size);
}

template <>
void Linear<API>::initSolver()
{
	_lin_solver.init(
		_K,
		_globalB,
		_localB,
		_lambda_map_sub_B1,
		_lambda_map_sub_B0,
		_lambda_map_sub_clst,
		_B1_duplicity,
		_f,
		_vec_c,
		_neighClusters
	);
}

}
