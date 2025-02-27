
#ifndef SRC_MATH_OPERATIONS_SORTING_PERMUTATION_H
#define SRC_MATH_OPERATIONS_SORTING_PERMUTATION_H

#include "math/primitives_new/vector_dense_view_new.h"
#include "math/primitives_new/permutation_view_new.h"



namespace espreso {
namespace math {
namespace operations {



template<typename T, typename I>
class sorting_permutation
{
public:
    void set_vector(VectorDenseView_new<T> * vec_);
    void set_permutation(PermutationView_new<I> * perm_);
    void perform();
    static void do_all(VectorDenseView_new<T> * vec, PermutationView_new<I> * perm);
private:
    VectorDenseView_new<T> * vec;
    PermutationView_new<I> * perm;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_SORTING_PERMUTATION_H */
