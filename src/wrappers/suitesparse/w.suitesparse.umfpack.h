
#ifndef SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_UMFPACK_H_
#define SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_UMFPACK_H_

#ifdef SUITESPARSE_HEADER_SUBDIR
#include "suitesparse/umfpack.h"
#else // SUITESPARSE_HEADER_DIRECT
#include "umfpack.h"
#endif

#include <complex>

namespace espreso {

template <typename I> inline void _defaults(double control[]);
template <> inline void _defaults<int>(double control[]) { umfpack_di_defaults(control); }
template <> inline void _defaults<long>(double control[]) { umfpack_dl_defaults(control); }

template <typename T, typename I> inline void _symbolic(const Matrix_CSC<T, I> &A, void **symbolic, double control[], double info[]);
template <> inline void _symbolic<double, int> (const Matrix_CSC<double, int>  &A, void **symbolic, double control[], double info[]) { umfpack_di_symbolic(A.nrows, A.nrows, A.cols, A.rows, nullptr, symbolic, control, info); }
template <> inline void _symbolic<double, long>(const Matrix_CSC<double, long> &A, void **symbolic, double control[], double info[]) { umfpack_dl_symbolic(A.nrows, A.nrows, A.cols, A.rows, nullptr, symbolic, control, info); }
template <> inline void _symbolic<std::complex<double>, int> (const Matrix_CSC<std::complex<double>, int>  &A, void **symbolic, double control[], double info[]) { umfpack_zi_symbolic(A.nrows, A.nrows, A.cols, A.rows, nullptr, nullptr, symbolic, control, info); }
template <> inline void _symbolic<std::complex<double>, long>(const Matrix_CSC<std::complex<double>, long> &A, void **symbolic, double control[], double info[]) { umfpack_zl_symbolic(A.nrows, A.nrows, A.cols, A.rows, nullptr, nullptr, symbolic, control, info); }

template <typename T, typename I> inline void _numeric(const Matrix_CSC<T, I> &A, void *symbolic, void **numeric, double control[], double info[]);
template <> inline void _numeric<double, int> (const Matrix_CSC<double, int>  &A, void *symbolic, void **numeric, double control[], double info[]) { umfpack_di_numeric(A.cols, A.rows, A.vals, symbolic, numeric, control, info); }
template <> inline void _numeric<double, long>(const Matrix_CSC<double, long> &A, void *symbolic, void **numeric, double control[], double info[]) { umfpack_dl_numeric(A.cols, A.rows, A.vals, symbolic, numeric, control, info); }
template <> inline void _numeric<std::complex<double>, int> (const Matrix_CSC<std::complex<double>, int>  &A, void *symbolic, void **numeric, double control[], double info[]) { umfpack_zi_numeric(A.cols, A.rows, reinterpret_cast<double*>(A.vals), nullptr, symbolic, numeric, control, info); }
template <> inline void _numeric<std::complex<double>, long>(const Matrix_CSC<std::complex<double>, long> &A, void *symbolic, void **numeric, double control[], double info[]) { umfpack_zl_numeric(A.cols, A.rows, reinterpret_cast<double*>(A.vals), nullptr, symbolic, numeric, control, info); }

template <typename T, typename I> inline void _solve(int sys, const Matrix_CSC<T, I> &A, T *x, const T *b, void *numeric, double control[], double info[]);
template <> inline void _solve<double, int> (int sys, const Matrix_CSC<double, int>  &A, double *x, const double *b, void *numeric, double control[], double info[]) { umfpack_di_solve(sys, A.cols, A.rows, A.vals, x, b, numeric, control, info); }
template <> inline void _solve<double, long>(int sys, const Matrix_CSC<double, long> &A, double *x, const double *b, void *numeric, double control[], double info[]) { umfpack_dl_solve(sys, A.cols, A.rows, A.vals, x, b, numeric, control, info); }
template <> inline void _solve<std::complex<double>, int> (int sys, const Matrix_CSC<std::complex<double>, int>  &A, std::complex<double> *x, const std::complex<double> *b, void *numeric, double control[], double info[]) { umfpack_zi_solve(sys, A.cols, A.rows, reinterpret_cast<double*>(A.vals), nullptr, reinterpret_cast<double*>(x), nullptr, reinterpret_cast<const double*>(b), nullptr, numeric, control, info); }
template <> inline void _solve<std::complex<double>, long>(int sys, const Matrix_CSC<std::complex<double>, long> &A, std::complex<double> *x, const std::complex<double> *b, void *numeric, double control[], double info[]) { umfpack_zl_solve(sys, A.cols, A.rows, reinterpret_cast<double*>(A.vals), nullptr, reinterpret_cast<double*>(x), nullptr, reinterpret_cast<const double*>(b), nullptr, numeric, control, info); }


}



#endif /* SRC_WRAPPERS_SUITESPARSE_W_SUITESPARSE_UMFPACK_H_ */
