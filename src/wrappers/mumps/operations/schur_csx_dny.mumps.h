
#ifndef SRC_WRAPPERS_MUMPS_OPERATIONS_SC_CSX_DNY_MUMPS_H_
#define SRC_WRAPPERS_MUMPS_OPERATIONS_SC_CSX_DNY_MUMPS_H_

#include "math/operations/schur_csx_dny.h"



namespace espreso {
namespace math {
namespace operations {



#ifdef HAVE_MUMPS



template<typename T, typename I>
struct schur_csx_dny_mumps_data;

template<typename T, typename I>
class schur_csx_dny_mumps : public schur_csx_dny<T,I>
{
public:
    schur_csx_dny_mumps();
    schur_csx_dny_mumps(const schur_csx_dny_mumps &) = delete;
    schur_csx_dny_mumps(schur_csx_dny_mumps &&) = default;
    schur_csx_dny_mumps & operator=(const schur_csx_dny_mumps &) = delete;
    schur_csx_dny_mumps & operator=(schur_csx_dny_mumps &&) = default;
    virtual ~schur_csx_dny_mumps();
public:
    virtual const char * get_name() override { return "schur_csx_dny_mumps"; }
protected:
    virtual void internal_preprocess() override;
    virtual void internal_perform_1() override;
    virtual void internal_perform_2() override;
    virtual void internal_solve_A11(VectorDenseView_new<T> & rhs, VectorDenseView_new<T> & sol) override;
private:
    std::unique_ptr<schur_csx_dny_mumps_data<T,I>> data;
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
class schur_csx_dny_mumps : public schur_csx_dny<T,I>
{
public:
    schur_csx_dny_mumps() { eslog::error("operation schur_csx_dny_mumps is not available\n"); }
    virtual ~schur_csx_dny_mumps() = default;
};



#endif



}
}
}

#endif /* SRC_WRAPPERS_MUMPS_OPERATIONS_SC_CSX_DNY_MUMPS_H_ */
