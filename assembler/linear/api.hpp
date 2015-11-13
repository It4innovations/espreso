
#include "linear.h"

namespace assembler {

template <>
void Linear<API>::KMf(size_t part, bool dynamics)
{
	_K[part] = this->_input.K;
}

template <>
void Linear<API>::RHS()
{
	_f[0] = std::vector<double>(this->_input.rhs.values, this->_input.rhs.values + this->_input.rhs.size);
	if (esconfig::MPIrank != 0) {
		return;
	}
	std::ofstream os("F_API.txt");
	eslocal index = 1;
	for (size_t i = 0; i < _f[0].size(); i++) {
		os << _f[0][i] << " ";
		if (index % 3 == 0) {
			os << std::endl;
		}
		index++;
	}
	os.close();
}

template <>
void Linear<API>::saveResult()
{
	std::stringstream ss;
	ss << "solution" << esconfig::MPIrank << ".txt";
	std::ofstream os(ss.str().c_str());
	eslocal index = 1;
	for (size_t i = 0; i < _prim_solution[0].size(); i++) {
		os << _prim_solution[0][i] << " ";
		if (index % 3 == 0) {
			os << std::endl;
		}
		index++;
	}
	os.close();
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
