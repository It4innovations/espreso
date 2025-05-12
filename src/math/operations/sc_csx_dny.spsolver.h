
#ifndef SRC_MATH_OPERATIONS_SC_CSX_DNY_SPSOLVER_H
#define SRC_MATH_OPERATIONS_SC_CSX_DNY_SPSOLVER_H

#include "math/operations/sc_csx_dny.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
struct sc_csx_dny_spsolver_data;

template<typename T, typename I>
class sc_csx_dny_spsolver : public sc_csx_dny<T,I>
{
public:
    sc_csx_dny_spsolver();
    sc_csx_dny_spsolver(const sc_csx_dny_spsolver &) = delete;
    sc_csx_dny_spsolver(sc_csx_dny_spsolver &&) = default;
    sc_csx_dny_spsolver & operator=(const sc_csx_dny_spsolver &) = delete;
    sc_csx_dny_spsolver & operator=(sc_csx_dny_spsolver &&) = default;
    virtual ~sc_csx_dny_spsolver();
protected:
    virtual void internal_preprocess() override;
    virtual void internal_perform_1() override;
    virtual void internal_perform_2() override;
    virtual void internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) override;
private:
    std::unique_ptr<sc_csx_dny_spsolver_data<T,I>> data;
protected:
    using sc_csx_dny<T,I>::A11;
    using sc_csx_dny<T,I>::A12;
    using sc_csx_dny<T,I>::A21;
    using sc_csx_dny<T,I>::A22;
    using sc_csx_dny<T,I>::A;
    using sc_csx_dny<T,I>::sc;
    using sc_csx_dny<T,I>::alpha;
    using sc_csx_dny<T,I>::need_solve_A11;
    using sc_csx_dny<T,I>::called_set_matrix;
    using sc_csx_dny<T,I>::called_preprocess;
    using sc_csx_dny<T,I>::called_perform;
    using sc_csx_dny<T,I>::size_matrix;
    using sc_csx_dny<T,I>::size_sc;
    using sc_csx_dny<T,I>::size_A11;
    using sc_csx_dny<T,I>::is_matrix_hermitian;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SC_CSX_DNY_SPSOLVER_H */
