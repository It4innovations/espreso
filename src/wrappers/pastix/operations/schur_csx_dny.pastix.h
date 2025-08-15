
#ifndef SRC_WRAPPERS_PASTIX_OPERATIONS_SC_CSX_DNY_PASTIX_H_
#define SRC_WRAPPERS_PASTIX_OPERATIONS_SC_CSX_DNY_PASTIX_H_

#include "math/operations/schur_csx_dny.h"



namespace espreso {
namespace math {
namespace operations {



#ifdef HAVE_PASTIX



template<typename T, typename I>
struct schur_csx_dny_pastix_data;

template<typename T, typename I>
class schur_csx_dny_pastix : public schur_csx_dny<T,I>
{
public:
    schur_csx_dny_pastix();
    schur_csx_dny_pastix(const schur_csx_dny_pastix &) = delete;
    schur_csx_dny_pastix(schur_csx_dny_pastix &&) = default;
    schur_csx_dny_pastix & operator=(const schur_csx_dny_pastix &) = delete;
    schur_csx_dny_pastix & operator=(schur_csx_dny_pastix &&) = default;
    virtual ~schur_csx_dny_pastix();
public:
    virtual const char * get_name() override { return "schur_csx_dny_pastix"; }
protected:
    virtual void internal_preprocess() override;
    virtual void internal_perform_1() override;
    virtual void internal_perform_2() override;
    virtual void internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) override;
private:
    std::unique_ptr<schur_csx_dny_pastix_data<T,I>> data;
protected:
    using schur_csx_dny<T,I>::A11;
    using schur_csx_dny<T,I>::A12;
    using schur_csx_dny<T,I>::A21;
    using schur_csx_dny<T,I>::A22;
    using schur_csx_dny<T,I>::A;
    using schur_csx_dny<T,I>::sc;
    using schur_csx_dny<T,I>::alpha;
    using schur_csx_dny<T,I>::need_solve_A11;
    using schur_csx_dny<T,I>::called_set_matrix;
    using schur_csx_dny<T,I>::called_preprocess;
    using schur_csx_dny<T,I>::called_perform;
    using schur_csx_dny<T,I>::size_matrix;
    using schur_csx_dny<T,I>::size_sc;
    using schur_csx_dny<T,I>::size_A11;
    using schur_csx_dny<T,I>::is_matrix_hermitian;
};



#else



template<typename T, typename I>
class schur_csx_dny_pastix : public schur_csx_dny<T,I>
{
public:
    schur_csx_dny_pastix() { eslog::error("operation schur_csx_dny_pastix is not available\n"); }
    virtual ~schur_csx_dny_pastix() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_PASTIX_OPERATIONS_SC_CSX_DNY_PASTIX_H_ */
