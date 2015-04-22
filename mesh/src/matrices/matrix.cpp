#include "matrix.h"

Matrix::Matrix(MatrixType type, MKL_INT rows, MKL_INT cols): _type(type), _rows(rows), _cols(cols)
{
	switch (_type) {
	case Matrix::SYMETRIC: {
		if (_rows != _cols) {
			fprintf(stderr, "Symmetric matrix has to have the same number of rows and columns.\n");
			exit(-1);
		}
	}
	}
}

void Matrix::transpose()
{
	if (_type == Matrix::SYMETRIC) {
		return;
	}

	makeTransposition();

	MKL_INT tmp = _rows;
	_rows = _cols;
	_cols = tmp;
}

std::ostream& operator<<(std::ostream& os, const Matrix &m)
{
	for(MKL_INT i = 0; i < m.rows(); i++) {
		for(MKL_INT j = 0; j < m.columns(); j++) {
			os << m(i, j) << " ";
		}
		os << std::endl;
	}
	os << std::endl;
	return os;
}
