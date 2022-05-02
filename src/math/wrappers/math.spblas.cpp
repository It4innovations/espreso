
#include "math/wrappers/math.spblas.h"
#include "math/wrappers/math.solver.h"
#include "esinfo/eslog.h"

#ifndef HAVE_MKL

namespace espreso {
namespace math {

template <>
void commit(Matrix_Dense<double> &x)
{

}

template <>
void commit(Matrix_Dense<std::complex<double> > &x)
{

}

template <>
void commit(Matrix_CSR<double> &x)
{

}

template <>
void commit(Matrix_CSR<std::complex<double> > &x)
{

}

template <>
void commit(Matrix_IJV<double> &x)
{

}

template <>
void commit(Matrix_IJV<std::complex<double> > &x)
{

}

template <>
void free(Matrix_Dense<double> &x)
{

}

template <>
void free(Matrix_Dense<std::complex<double> > &x)
{

}

template <>
void free(Matrix_CSR<double> &x)
{

}

template <>
void free(Matrix_CSR<std::complex<double> > &x)
{

}

template <>
void free(Matrix_IJV<double> &x)
{

}

template <>
void free(Matrix_IJV<std::complex<double> > &x)
{

}


template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_CSR<double> &a, const double &beta, const Vector_Dense<double> &x)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_CSR<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

template <>
void apply(Vector_Dense<double> &y, const double &alpha, const Matrix_IJV<double> &a, const double &beta, const Vector_Dense<double> &x)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

template <>
void apply(Vector_Dense<std::complex<double> > &y, const std::complex<double> &alpha, const Matrix_IJV<std::complex<double> > &a, const std::complex<double> &beta, const Vector_Dense<std::complex<double> > &x)
{
	eslog::error("calling of empty SpBLAS wrapper.\n");
}

}
}

#endif
