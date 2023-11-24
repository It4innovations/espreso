
#if defined(MY_HIP) || defined(MY_ROC)

#include <cstdio>

#include <hip/hip_runtime.h>

#include "my_common.hpp"
#include "matrices.hpp"



// for compatibility with older versions (< rocm-5.4)
#ifdef MY_OLD_HIP
#define HIPSPARSE_ORDER_COL HIPSPARSE_ORDER_COLUMN
#endif





template<typename T>
class my_hipd_allocator
{
public:
    static constexpr bool is_data_host_accessible = false;
public:
    using value_type = T;
    int device;
    my_hipd_allocator()
    {
        CHECK(hipGetDevice(&device));
    }
    explicit my_hipd_allocator(int d) : device(d) { }
    T * allocate(size_t count)
    {
        T * ptr;
        int orig_device;
        CHECK(hipGetDevice(&orig_device));
        CHECK(hipSetDevice(device));
        CHECK(hipMalloc(&ptr, count * sizeof(T)));
        CHECK(hipSetDevice(orig_device));
        return ptr;
    }
    void deallocate(T * ptr, size_t count)
    {
        int orig_device;
        CHECK(hipGetDevice(&orig_device));
        CHECK(hipSetDevice(device));
        CHECK(hipFree(ptr));
        CHECK(hipSetDevice(orig_device));
    }
};

template<typename T, typename I>
using MatrixDense_hipd = MatrixDense<T,I,my_hipd_allocator>;
template<typename T, typename I>
using MatrixCSR_hipd = MatrixCSR<T,I,my_hipd_allocator>;
template<typename I>
using Permutation_hipd = Permutation<I,my_hipd_allocator>;



template<typename T>
class my_hiph_allocator
{
public:
    static constexpr bool is_data_host_accessible = true;
public:
    using value_type = T;
    my_hiph_allocator() = default;
    T * allocate(size_t count)
    {
        T * ptr;
        CHECK(hipHostMalloc(&ptr, count * sizeof(T)));
        return ptr;
    }
    void deallocate(T * ptr, size_t count)
    {
        CHECK(hipHostFree(ptr));
    }
};

template<typename T, typename I>
using MatrixDense_hiph = MatrixDense<T,I,my_hiph_allocator>;
template<typename T, typename I>
using MatrixCSR_hiph = MatrixCSR<T,I,my_hiph_allocator>;
template<typename I>
using Permutation_hiph = Permutation<I,my_hiph_allocator>;

template<typename T, typename I>
using MatrixDense_l_hiph = MatrixDense_l<T,I,my_hiph_allocator>;
template<typename T, typename I>
using MatrixCSR_l_hiph = MatrixCSR_l<T,I,my_hiph_allocator>;
template<typename I>
using Permutation_l_hiph = Permutation_l<I,my_hiph_allocator>;



template<typename T>
class my_hipm_allocator
{
public:
    static constexpr bool is_data_host_accessible = true;
public:
    using value_type = T;
    my_hipm_allocator() = default;
    T * allocate(size_t count)
    {
        T * ptr;
        CHECK(hipMallocManaged(&ptr, count * sizeof(T)));
        return ptr;
    }
    void deallocate(T * ptr, size_t count)
    {
        CHECK(hipFree(ptr));
    }
};

template<typename T, typename I>
using MatrixDense_hipm = MatrixDense<T,I,my_hipm_allocator>;
template<typename T, typename I>
using MatrixCSR_hipm = MatrixCSR<T,I,my_hipm_allocator>;
template<typename I>
using Permutation_hipm = Permutation<I,my_hipm_allocator>;





template<typename T, typename I, template<typename> typename A1, template<typename> typename A2>
void copy_matrix_submit(MatrixDense<T,I,A1> & output, const MatrixDense<T,I,A2> & input, hipMemcpyKind direction, hipStream_t stream = 0)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols) MY_ABORT("copy_matrix_submit: output matrix has wrong dimensions");
    CHECK(hipMemcpy2DAsync(output.vals, output.ld * sizeof(T), input.vals, input.ld * sizeof(T), input.ncols * sizeof(T), input.nrows, direction, stream));
}

template<typename T, typename I>
void copy_matrix_submit(MatrixDense_hiph<T,I> & output, const MatrixDense_hipd<T,I> & input, hipStream_t stream = 0)
{
    copy_matrix_submit(output, input, hipMemcpyDeviceToHost, stream);
}

template<typename T, typename I>
void copy_matrix_submit(MatrixDense_hipd<T,I> & output, const MatrixDense_hiph<T,I> & input, hipStream_t stream = 0)
{
    copy_matrix_submit(output, input, hipMemcpyHostToDevice, stream);
}

template<typename T, typename I>
void copy_matrix_submit(MatrixDense_hipd<T,I> & output, const MatrixDense_hipd<T,I> & input, hipStream_t stream = 0)
{
    copy_matrix_submit(output, input, hipMemcpyDeviceToDevice, stream);
}

template<typename T, typename I>
void copy_matrix_d2h(MatrixDense<T,I> & output, const MatrixDense_hipd<T,I> & input)
{
    copy_matrix_submit(output, input, hipMemcpyDeviceToHost);
    CHECK(hipDeviceSynchronize());
}



template<typename T, typename I, template<typename> typename A1, template<typename> typename A2>
void copy_matrix_submit(MatrixCSR<T,I,A1> & output, const MatrixCSR<T,I,A2> & input, hipMemcpyKind direction, hipStream_t stream = 0, bool copy_rowptrs = true, bool copy_colidxs = true, bool copy_vals = true)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != input.nvals) MY_ABORT("copy_matrix_submit: output matrix has wrong dimensions");
    if(copy_rowptrs) CHECK(hipMemcpyAsync(output.rowptrs, input.rowptrs, (input.nrows+1) * sizeof(I), direction, stream));
    if(copy_colidxs) CHECK(hipMemcpyAsync(output.colidxs, input.colidxs,  input.nvals    * sizeof(I), direction, stream));
    if(copy_vals)    CHECK(hipMemcpyAsync(output.vals,    input.vals,     input.nvals    * sizeof(T), direction, stream));
}

template<typename T, typename I>
void copy_matrix_submit(MatrixCSR_hiph<T,I> & output, const MatrixCSR_hipd<T,I> & input, hipStream_t stream = 0, bool copy_rowptrs = true, bool copy_colidxs = true, bool copy_vals = true)
{
    copy_matrix_submit(output, input, hipMemcpyDeviceToHost, stream, copy_rowptrs, copy_colidxs, copy_vals);
}

template<typename T, typename I>
void copy_matrix_submit(MatrixCSR_hipd<T,I> & output, const MatrixCSR_hiph<T,I> & input, hipStream_t stream = 0, bool copy_rowptrs = true, bool copy_colidxs = true, bool copy_vals = true)
{
    copy_matrix_submit(output, input, hipMemcpyHostToDevice, stream, copy_rowptrs, copy_colidxs, copy_vals);
}

template<typename T, typename I>
void copy_matrix_submit(MatrixCSR_hipd<T,I> & output, const MatrixCSR_hipd<T,I> & input, hipStream_t stream = 0, bool copy_rowptrs = true, bool copy_colidxs = true, bool copy_vals = true)
{
    copy_matrix_submit(output, input, hipMemcpyDeviceToDevice, stream, copy_rowptrs, copy_colidxs, copy_vals);
}



template<typename I, template<typename> typename A1, template<typename> typename A2>
void copy_perm_submit(Permutation<I,A1> & output, const Permutation<I,A2> & input, hipMemcpyKind direction, hipStream_t stream = 0)
{
    if(output.size != input.size) MY_ABORT("copy_perm_submit: output matrix has wrong dimensions");
    CHECK(hipMemcpyAsync(output.forward,  input.forward,  input.size * sizeof(I), direction, stream));
    CHECK(hipMemcpyAsync(output.backward, input.backward, input.size * sizeof(I), direction, stream));
}

template<typename I>
void copy_perm_submit(Permutation_hiph<I> & output, const Permutation_hipd<I> & input, hipStream_t stream = 0)
{
    copy_perm_submit(output, input, hipMemcpyDeviceToHost, stream);
}

template<typename I>
void copy_perm_submit(Permutation_hipd<I> & output, const Permutation_hiph<I> & input, hipStream_t stream = 0)
{
    copy_perm_submit(output, input, hipMemcpyHostToDevice, stream);
}

template<typename I>
void copy_perm_submit(Permutation_hipd<I> & output, const Permutation_hipd<I> & input, hipStream_t stream = 0)
{
    copy_perm_submit(output, input, hipMemcpyDeviceToDevice, stream);
}





template<typename T, typename I>
__global__ void my_submatrix_kernel(T * densevals, I ld, const I * rowptrs, const I * colidxs, const T * sparsevals, I rowoffset, I colstart, I colend)
{
    // gridDim.x = number of rows of output matrix
    I r_out = blockIdx.x;
    I r_in = blockIdx.x + rowoffset;
    I start = rowptrs[r_in];
    I end = rowptrs[r_in+1];
    T * rowvals = densevals + r_out * ld;

    for(I i = start + threadIdx.x; i < end; i += blockDim.x)
    {
        I c_in = colidxs[i];
        I c_out = c_in - colstart;
        if(c_in >= colend) break;
        if(c_in >= colstart) rowvals[c_out] += sparsevals[i];
    }
}

template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
void submit_submatrix(MatrixDense<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, I rowstart, I rowend, I colstart, I colend, hipStream_t stream)
{
    if(output.nrows != rowend - rowstart || output.ncols != colend - colstart) MY_ABORT("submit_submatrix: output matrix has wrong dimensions");

    CHECK(hipMemsetAsync(output.vals, 0, output.size_alloc() * sizeof(T), stream));

    I nrows_out = rowend - rowstart;
    int tpb = 32;
    int bpg = nrows_out;
    my_submatrix_kernel<<< bpg, tpb, 0, stream >>>(output.vals, output.ld, input.rowptrs, input.colidxs, input.vals, rowstart, colstart, colend);
    CHECK(hipPeekAtLastError());
}





template<typename T, typename I>
__global__ void my_permute_matrix_rows_kernel(T * output, I out_ld, const T * input, I in_ld, const I * perm, I nrows, I ncols)
{
    // made primarily for 1-column matrices
    I row_in = blockIdx.x * blockDim.x + threadIdx.x;
    if(row_in >= nrows) return;
    I row_out = perm[row_in];

    const T * rin = input + row_in * in_ld;
    T * rout = output + row_out * out_ld;

    for(I c = 0; c < ncols; c++)
    {
        rout[c] = rin[c];
    }
}

template<typename T, typename I>
void submit_permute_matrix_rows(MatrixDense_hipd<T,I> & output, const MatrixDense_hipd<T,I> & input, const Permutation_hipd<I> & perm, bool invperm = false, hipStream_t stream = 0)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols) MY_ABORT("submit_permute_matrix_rows: output matrix has wrong dimensions");

    const I * permvals;
    if(invperm) permvals = perm.backward;
    else        permvals = perm.forward;
    int tpb = 32;
    int bpg = (input.nrows - 1) / tpb + 1;
    my_permute_matrix_rows_kernel<<< bpg, tpb, 0, stream >>>(output.vals, output.ld, input.vals, input.ld, permvals, input.nrows, input.ncols);
    CHECK(hipPeekAtLastError());
}





__global__ void dummy_init_kernel(int x)
{
    if(x == 17) printf("it is 17. %d %d\n", (int)blockIdx.x, (int)threadIdx.x);
}

void run_dummy_init_kernel(int x)
{
    dummy_init_kernel<<<1,1>>>(x);
    CHECK(hipDeviceSynchronize());
}





template<typename T, typename I>
static void print_matrix(const MatrixDense_hipd<T,I> & Md, const char * name)
{
    MatrixDense_hiph<T,I> Mh(Md.nrows, Md.ncols, Md.ld, true);
    copy_matrix_submit(Mh, Md);
    CHECK(hipDeviceSynchronize());
    print_matrix(Mh, name);
}

template<typename T, typename I>
static void print_matrix(const MatrixCSR_hipd<T,I> & Md, const char * name)
{
    MatrixCSR_hiph<T,I> Mh(Md.nrows, Md.ncols, Md.nvals, true);
    copy_matrix_submit(Mh, Md);
    CHECK(hipDeviceSynchronize());
    print_matrix(Mh, name);
}





template<typename C>
void my_hip_submit_host_function(hipStream_t stream, C && c)
{
    C * cp = new C(std::move(c));

#ifdef MY_OLD_HIP
    CHECK(hipStreamAddCallback(stream, [](hipStream_t stream_, hipError_t status_, void * arg){
        C * cpi = reinterpret_cast<C*>(arg);
        (*cpi)();
        delete cpi;
    }, cp, 0));
#else
    CHECK(hipLaunchHostFunc(stream, [](void * arg){
        C * cpi = reinterpret_cast<C*>(arg);
        (*cpi)();
        delete cpi;
    }, cp));
#endif
}



void hip_malloc_max_memory(char ** memory, size_t * memory_size_B)
{
    size_t starting_size;
    double coef_percent = 95;

    {
        int device;
        CHECK(hipGetDevice(&device));
        hipDeviceProp_t props;
        CHECK(hipGetDeviceProperties(&props, device));
        starting_size = props.totalGlobalMem;
    }

    for(size_t memsize = starting_size; memsize > 0; memsize = memsize * coef_percent / 100)
    {
        char * mem;
        hipError_t err = hipMalloc(&mem, memsize);

        if(err == hipSuccess)
        {
            *memory = mem;
            *memory_size_B = memsize;
            return;
        }
        if(err != hipErrorMemoryAllocation) CHECK(err);
    }

    throw std::runtime_error("could not allocate any memory");
}



#endif
