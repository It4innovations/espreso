
#ifndef SRC_MATH_OPERATIONS_AUXILIARY_TRI_PARTITION_HERK_H
#define SRC_MATH_OPERATIONS_AUXILIARY_TRI_PARTITION_HERK_H



class tri_partition_herk
{
public:
    void set_config(char algorithm_, char direction_, int parameter_, char herk_strategy_);
    void set_system(size_t size_n_, size_t size_k_);
    void setup();
    size_t get_num_chunks();
    void set_output_partition(VectorDenseView_new<size_t> * partition_);
    void perform();
private:
    char algorithm = '_'; // Uniform, Minimal work
    char direction = '_'; // along N direction, along K direction
    char herk_strategy = '_'; // sTairs, sQuares
    int parameter = 0;
    size_t size_n = 0;
    size_t size_k = 0;
    size_t partition_range = 0;
    size_t num_chunks = 0;
    VectorDenseView_new<size_t> * partition = nullptr;
    bool set_config_called = false;
    bool set_system_called = false;
    bool setup_called = false;
};



#endif /* SRC_MATH_OPERATIONS_AUXILIARY_TRI_PARTITION_HERK_H */
