#include "matrix.h"

std::ostream& operator<<(std::ostream& os, const Matrix &m)
{
	for (size_t i = 0; i < m.rows(); i++) {
		for (size_t j = 0; j < m.columns(); j++) {
			os << m(i, j) << " ";
		}
		os << std::endl;
	}
	os << std::endl;
	return os;
}
