
#ifdef HAVE_MKL

#include "esinfo/eslog.h"

namespace espreso {

template <typename T, typename I>
inline void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T,I> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
{
	eslog::error("MKL PARDISO does not provide factors.\n");
}

template <typename T, typename I>
inline void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T,I> &/*U*/, bool /*copyPattern*/, bool /*copyValues*/)
{
	eslog::error("MKL PARDISO does not provide factors.\n");
}

}

#endif
