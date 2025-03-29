
#ifndef SRC_MATH_OPERATIONS_AUXILIARY_TRI_PARTITION_TRSM_H
#define SRC_MATH_OPERATIONS_AUXILIARY_TRI_PARTITION_TRSM_H

#include "math/primitives_new/vector_dense_view_new.h"



namespace espreso {
namespace math {
namespace operations {



class tri_partition_trsm
{
public:
    tri_partition_trsm() = default;
    tri_partition_trsm(const tri_partition_trsm &) = delete;
    tri_partition_trsm(tri_partition_trsm &&) = default;
    tri_partition_trsm & operator=(const tri_partition_trsm &) = delete;
    tri_partition_trsm & operator=(tri_partition_trsm &&) = default;
    ~tri_partition_trsm() = default;
public:
    void set_config(char algorithm_, char direction_, int parameter_);
    void set_system(size_t sys_size_, size_t sys_nrhs_);
    void setup();
    size_t get_num_chunks();
    void set_output_partition(VectorDenseView_new<size_t> * partition_);
    void perform();
private:
    char algorithm = '_'; // Uniform, Minimal work
    char direction = '_'; // Horizontal, Vertical
    int parameter = 0;
    size_t sys_size = 0;
    size_t sys_nrhs = 0;
    size_t partition_range = 0;
    size_t num_chunks = 0;
    VectorDenseView_new<size_t> * partition = nullptr;
    bool set_config_called = false;
    bool set_system_called = false;
    bool setup_called = false;
};



}
}
}

#endif /* SRC_MATH_OPERATIONS_AUXILIARY_TRI_PARTITION_TRSM_H */
