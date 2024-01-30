
#ifdef HAVE_SUITESPARSE

#include <algorithm>
#include "w.suitesparse.cholmod.h"
#include "esinfo/eslog.h"

namespace espreso {

template<typename T, typename I>
struct Solver_External_Representation {
	cholmod_common cm_common;
	cholmod_factor * cm_factor_super = nullptr;
	cholmod_factor * cm_factor_simpl = nullptr;
	cholmod_sparse * cm_matrix_view = nullptr;
    const Matrix_CSR<T, I> * matrix = nullptr;
	Vector_Dense<I> map_simpl_super;
	char zerodrop;
	int stage = 0; // 0 = completely uninitialized, 1 = initialized without matrix, 2 = have matrix, 3 = symbolic factorization done, 4 = numeric factorization done
};



template <typename T, typename I>
template <typename A>
inline void DirectSparseSolver<T, I>::getFactorL(Matrix_CSR<T,I,A> &/*L*/, bool /*copyPattern*/, bool /*copyValues*/)
{
	eslog::error("L factor is not provided\n");
}

template <typename T, typename I>
template <typename A>
inline void DirectSparseSolver<T, I>::getFactorU(Matrix_CSR<T,I,A> &U, bool copyPattern, bool copyValues)
{
	if(ext->stage < 3) eslog::error("getFactorU: invalid order of operations in spsolver\n");
    if(copyValues && ext->stage < 4) eslog::error("getFactorU: invalid order of operations in spsolver\n");
    if((size_t)U.nrows != ext->cm_factor_simpl->n || (size_t)U.ncols != ext->cm_factor_simpl->n) eslog::error("getFactorU: output matrix has wrong dimensions\n");

	U.resize(ext->cm_factor_simpl->n, ext->cm_factor_simpl->n, ext->cm_factor_simpl->nzmax);

    if(copyPattern) std::copy_n(static_cast<I*>(ext->cm_factor_simpl->p), ext->cm_factor_simpl->n+1, U.rows);
    if(copyPattern) std::copy_n(static_cast<I*>(ext->cm_factor_simpl->i), ext->cm_factor_simpl->nzmax, U.cols);
    if(copyValues) for(I i = 0; i < ext->map_simpl_super.size; i++) U.vals[i] = reinterpret_cast<T*>(ext->cm_factor_super->x)[ext->map_simpl_super.vals[i]];
}

}

#endif
