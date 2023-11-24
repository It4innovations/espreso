
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <limits>
#include <memory>
#include <deque>
#include <stdexcept>

#include "my_timer.hpp"

template<typename T, typename I, template<typename> typename A = my_stdallocator_wrapper>
class MatrixCSR
{
public:
    using float_type = T;
    using int_type = I;
    using allocator_type_T = A<T>;
    using allocator_type_I = A<I>;
public:
    T * vals = nullptr;
    I * colidxs = nullptr;
    I * rowptrs = nullptr;
    A<T> at;
    A<I> ai;
    I nrows;
    I ncols;
    I nvals;
    bool should_i_free;
public:
    MatrixCSR() : MatrixCSR(0, 0, 0, false) { }
    MatrixCSR(const A<T> & at_, const A<I> & ai_) : MatrixCSR(0, 0, 0, false, at_, ai_) { }
    MatrixCSR(I nrows_, I ncols_, I nvals_, bool do_allocate_, const A<T> & at_ = A<T>(), const A<I> & ai_ = A<I>())
        : nrows{nrows_}, ncols{ncols_}, nvals{nvals_}, should_i_free{do_allocate_}, at{at_}, ai{ai_}
    {
        if(do_allocate_)
        {
            allocate();
        }
    }
    MatrixCSR(const MatrixCSR & other) = delete;
    MatrixCSR(MatrixCSR && other)
        : at(other.at), ai(other.ai), nrows(other.nrows), ncols(other.ncols), nvals(other.nvals), vals(other.vals), colidxs(other.colidxs), rowptrs(other.rowptrs), should_i_free(other.should_i_free)
    {
        other.should_i_free = false;
    }
    ~MatrixCSR()
    {
        deallocate();
    }
    MatrixCSR & operator=(const MatrixCSR & other) = delete;
    MatrixCSR & operator=(MatrixCSR && other)
    {
        deallocate();
        at = other.at;
        ai = other.ai;
        nrows = other.nrows;
        ncols = other.ncols;
        nvals = other.nvals;
        vals = other.vals;
        colidxs = other.colidxs;
        rowptrs = other.rowptrs;
        should_i_free = other.should_i_free;
        if(&other != this) other.should_i_free = false;
        return *this;
    }
    void resize(I nrows_, I ncols_, I nvals_, bool do_allocate_)
    {
        nrows = nrows_;
        ncols = ncols_;
        nvals = nvals_;
        should_i_free = do_allocate_;
        if(do_allocate_)
        {
            allocate();
        }
    }
    void allocate()
    {
        deallocate();
        if(nvals > 0) vals = at.allocate(nvals);
        if(nvals > 0) colidxs = ai.allocate(nvals);
        rowptrs = ai.allocate(nrows + 1);
        should_i_free = true;
    }
    void deallocate()
    {
        if(should_i_free)
        {
            at.deallocate(vals, nvals);
            ai.deallocate(colidxs, nvals);
            ai.deallocate(rowptrs, nrows + 1);
            vals = nullptr;
            colidxs = nullptr;
            rowptrs = nullptr;
        }
        should_i_free = false;
    }
};



template<typename T, typename I, template<typename> typename A = my_stdallocator_wrapper>
class MatrixDense
{
public:
    using float_type = T;
    using int_type = I;
    using allocator_type_T = A<T>;
    using allocator_type_I = A<I>;
public:
    T * vals = nullptr;
    A<T> at;
    I nrows;
    I ncols;
    I ld;
    bool should_i_free;
public:
    MatrixDense() : MatrixDense(0, 0, -1, false) { }
    MatrixDense(const A<T> & at_) : MatrixDense(0, 0, -1, false, at_) { }
    MatrixDense(I nrows_, I ncols_, I ld_, bool do_allocate_, const A<T> & at_ = A<T>())
        : nrows{nrows_}, ncols{ncols_}, ld{ld_}, should_i_free{do_allocate_}, at{at_}
    {
        if(ld < 0) ld = ncols;
        if(do_allocate_) allocate();
    }
    MatrixDense(const MatrixDense & other) = delete;
    MatrixDense(MatrixDense && other)
        : at(other.at), vals(other.vals), nrows(other.nrows), ncols(other.ncols), ld(other.ld), should_i_free(other.should_i_free)
    {
        other.should_i_free = false;
    }
    ~MatrixDense()
    {
        deallocate();
    }
    MatrixDense & operator=(const MatrixDense & other) = delete;
    MatrixDense & operator=(MatrixDense && other)
    {
        deallocate();
        at = other.at;
        vals = other.vals;
        nrows = other.nrows;
        ncols = other.ncols;
        ld = other.ld;
        should_i_free = other.should_i_free;
        if(&other != this) other.should_i_free = false;
        return *this;
    }
    void resize(I nrows_, I ncols_, I ld_, bool do_allocate_)
    {
        if(ld_ < 0) ld_ = ncols_;

        nrows = nrows_;
        ncols = ncols_;
        ld = ld_;
        should_i_free = do_allocate_;
        if(do_allocate_) allocate();
    }
    void allocate()
    {
        deallocate();
        if(size_alloc() != ((size_t)nrows * ld)) MY_ABORT("MatrixDense::allocate: matrix too large for the used index_type");
        if(size_alloc() > 0) vals = at.allocate(size_alloc());
        should_i_free = true;
    }
    void deallocate()
    {
        if(should_i_free)
        {
            at.deallocate(vals, size_alloc());
            vals = nullptr;
        }
        should_i_free = false;
    }
    I size_alloc() const
    {
        return nrows * ld;
    }
    I size_matrix() const
    {
        return nrows * ncols;
    }
};



enum struct PermutationDirection
{
    Forward,
    Backward
};



template<typename I, template<typename> typename A = my_stdallocator_wrapper>
class Permutation
{
public:
    using int_type = I;
    using allocator_type_I = A<I>;
public:
    I * forward = nullptr;
    I * backward = nullptr;
    A<I> ai;
    I size;
    bool should_i_free;
public:
    Permutation() : Permutation(0, false) { }
    Permutation(const A<I> & ai_) : Permutation(0, false, ai_) { }
    Permutation(I size_, bool do_allocate_, const A<I> & ai_ =  A<I>())
        : size{size_}, should_i_free{do_allocate_}, ai{ai_}
    {
        if(do_allocate_)
        {
            allocate();
        }
    }
    Permutation(const Permutation & other) = delete;
    Permutation(Permutation && other)
        : ai(other.ai), forward(other.forward), backward(other.backward), size(other.size), should_i_free(other.should_i_free)
    {
        other.should_i_free = false;
    }
    ~Permutation()
    {
        deallocate();
    }
    Permutation & operator=(const Permutation & other) = delete;
    Permutation & operator=(Permutation && other)
    {
        deallocate();
        ai = other.ai;
        forward = other.forward;
        backward = other.backward;
        size = other.size;
        should_i_free = other.should_i_free;
        if(&other != this) other.should_i_free = false;
        return *this;
    }
    void resize(I size_, bool do_allocate_)
    {
        size = size_;
        should_i_free = do_allocate_;
        if(do_allocate_)
        {
            allocate();
        }
    }
    void allocate()
    {
        deallocate();
        if(size > 0) forward = ai.allocate(size);
        if(size > 0) backward = ai.allocate(size);
        should_i_free = true;
    }
    void deallocate()
    {
        if(should_i_free)
        {
            ai.deallocate(forward, size);
            ai.deallocate(backward, size);
            forward = nullptr;
            backward = nullptr;
        }
        should_i_free = false;
    }
    void invert()
    {
        std::swap(forward, backward);
    }
};

template<typename T, typename I, template<typename> typename A = my_stdallocator_wrapper>
using MatrixCSR_l = MatrixCSR<T,I,my_limited_allocator_outer<A>::template inner>;

template<typename T, typename I, template<typename> typename A = my_stdallocator_wrapper>
using MatrixDense_l = MatrixDense<T,I,my_limited_allocator_outer<A>::template inner>;

template<typename I, template<typename> typename A = my_stdallocator_wrapper>
using Permutation_l = Permutation<I,my_limited_allocator_outer<A>::template inner>;



struct putimers
{
    my_timer read, sort, write;
};

struct permtimers
{
    my_timer upper2full, vecinit, traverse, invperm, doperm;
    putimers pu;
};

struct sctimers
{
    my_timer perm, dosc, copy;
    my_timer subS, boundcalc, subX, mkl, trf, trs, syrk;
    permtimers pt;
};






























template<typename T, typename I, template<typename> typename A>
static bool save_matrix(const MatrixDense<T,I,A> & M, const char * filepath)
{
    FILE * f = fopen(filepath, "w");
    if(f == nullptr)
    {
        fprintf(stderr, "Could not open file '%s'\n", filepath);
        return false;
    }

    fprintf(f, "%lld %lld\n", (long long)M.nrows, (long long)M.ncols);
    for(I r = 0; r < M.nrows; r++)
    {
        for(I c = 0; c < M.ncols; c++)
        {
            fprintf(f, "%+.15e ", (double)M.vals[r * M.ncols + c]);
        }
        fprintf(f, "\n");
    }

    fclose(f);

    return true;
}



template<typename T, typename I, template<typename> typename A>
static bool save_matrix(const MatrixCSR<T,I,A> & M, const char * filepath)
{
    FILE * f = fopen(filepath, "w");
    if(f == nullptr)
    {
        fprintf(stderr, "Could not open file '%s'\n", filepath);
        return false;
    }

    fprintf(f, "%lld %lld %lld\n", (long long)M.nrows, (long long)M.ncols, (long long)M.nvals);
    for(I r = 0; r < M.nrows; r++)
    {
        I start = M.rowptrs[r];
        I end = M.rowptrs[r+1];
        for(I i = start; i < end; i++)
        {
            I c = M.colidxs[i];
            double v = static_cast<double>(M.vals[i]);
            fprintf(f, "%6lld %6lld %+25.15e\n", (long long)r, (long long)c, (double)v);
        }
    }

    fclose(f);

    return true;
}



template<typename T, typename I, template<typename> typename A>
static bool load_matrix_csr(MatrixCSR<T,I,A> * output, const char * filepath)
{
    FILE * f = fopen(filepath, "r");
    if(f == nullptr)
    {
        fprintf(stderr, "Could not open file '%s'\n", filepath);
        return false;
    }

    I nrows, ncols, nvals;
    fscanf(f, std::is_same_v<I,int32_t> ? "%d%d%d" : "%ld%ld%ld", &nrows, &ncols, &nvals);
    output->resize(nrows, ncols, nvals, true);
    I lastrow = 0;
    output->rowptrs[0] = 0;
    for(I i = 0; i < output->nvals; i++)
    {
        I row, col;
        T val;
        fscanf(f, std::is_same_v<I,int32_t> ? "%d%d" : "%ld%ld", &row, &col);
        fscanf(f, std::is_same_v<T,float> ? "%f" : "%lf",  &val);
        while(row > lastrow)
        {
            lastrow++;
            output->rowptrs[lastrow] = i;
        }
        output->vals[i] = val;
        output->colidxs[i] = col;
    }
    while(nrows > lastrow)
    {
        lastrow++;
        output->rowptrs[lastrow] = nvals;
    }
    output->rowptrs[output->nrows] = output->nvals;

    fclose(f);

    return true;
}



template<typename T, typename I, template<typename> typename A>
static bool load_matrix_dense(MatrixDense<T,I,A> * output, const char * filepath)
{
    FILE * f = fopen(filepath, "r");
    if(f == nullptr)
    {
        fprintf(stderr, "Could not open file '%s'\n", filepath);
        return false;
    }

    I nrows, ncols;
    fscanf(f, std::is_same_v<I,int32_t> ? "%d%d" : "%ld%ld", &nrows, &ncols);
    output->resize(nrows, ncols, -1, true);
    for(I r = 0; r < nrows; r++)
    {
        for(I c = 0; c < ncols; c++)
        {
            fscanf(f, std::is_same_v<T,float> ? "%f" : "%lf", output->vals + r * output->ld + c);
        }
    }

    fclose(f);

    return true;
}



template<typename T, typename I, template<typename> typename A>
static bool load_vector_dense(MatrixDense<T,I,A> * output, const char * filepath)
{
    FILE * f = fopen(filepath, "r");
    if(f == nullptr)
    {
        fprintf(stderr, "Could not open file '%s'\n", filepath);
        return false;
    }

    std::vector<T> vals;
    while(true)
    {
        int nscanned;
        T currval;
        if constexpr(std::is_same_v<T,int>)       nscanned = fscanf(f, "%d",   &currval);
        if constexpr(std::is_same_v<T,long>)      nscanned = fscanf(f, "%ld",  &currval);
        if constexpr(std::is_same_v<T,long long>) nscanned = fscanf(f, "%lld", &currval);
        if constexpr(std::is_same_v<T,float>)     nscanned = fscanf(f, "%f",   &currval);
        if constexpr(std::is_same_v<T,double>)    nscanned = fscanf(f, "%f",   &currval);
        if(nscanned <= 0)
            break;
        vals.push_back(currval);
    }

    fclose(f);

    output->resize(vals.size(), 1, true);
    std::copy(vals.begin(), vals.end(), output->vals);

    return true;
}



template<typename T, typename I, template<typename> typename A>
static bool contains_nan(MatrixDense<T,I,A> & M)
{
    static_assert(A<T>::is_data_host_accessible, "Matrix data has to be host accessible");
    for(size_t i = 0; i < M.size_alloc(); i++) if(my_isnan(M.vals[i])) return true;
    return false;
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void sparse_to_dense(MatrixDense<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols) MY_ABORT("sparse_to_dense: output matrix has wrong dimensions");

    std::fill(output.vals, output.vals + output.size_alloc(), T{0});
    for(I r = 0; r < input.nrows; r++)
    {
        I i_start = input.rowptrs[r];
        I i_end = input.rowptrs[r+1];
        for(I i = i_start; i < i_end; i++)
        {
            I c = input.colidxs[i];
            output.vals[r * output.ld + c] = input.vals[i];
        }
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void sparse_to_dense_transpose(MatrixDense<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input)
{
    if(output.nrows != input.ncols || output.ncols != input.nrows) MY_ABORT("sparse_to_dense_transpose: output matrix has wrong dimensions");

    std::fill(output.vals, output.vals + output.size_alloc(), T{0});
    for(I r = 0; r < input.nrows; r++)
    {
        I i_start = input.rowptrs[r];
        I i_end = input.rowptrs[r+1];
        for(I i = i_start; i < i_end; i++)
        {
            I c = input.colidxs[i];
            output.vals[c * output.ld + r] = input.vals[i];
        }
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai, template<typename> typename Ap>
static void sparse_to_dense_permcol(MatrixDense<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, const Permutation<I,Ap> & perm)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols) MY_ABORT("sparse_to_dense_permcol: output matrix has wrong dimensions");

    std::fill(output.vals, output.vals + output.nrows * output.ncols, 0);

    for(I r = 0; r < input.nrows; r++)
    {
        I row = r;
        I i_start = input.rowptrs[r];
        I i_end = input.rowptrs[r+1];
        for(I i = i_start; i < i_end; i++)
        {
            I col = perm.backward[input.colidxs[i]];
            output.vals[row * output.ncols + col] = input.vals[i];
        }
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void sparse_to_dense_select_rows(MatrixDense<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, const std::vector<I> & rows)
{
    if(output.nrows != rows.size() || output.ncols != input.ncols) MY_ABORT("sparse_to_dense_select_rows: output matrix has wrong dimensions");

    std::fill(output.vals, output.vals + output.get_nvals(), I{0});

    for(I r_out = 0; r_out < output.nrows; r_out++)
    {
        I r_in = rows[r_out];
        I start = input.rowptrs[r_in];
        I end = input.rowptrs[r_in+1];
        for(I i = start; i < end; i++)
        {
            I c = input.colidxs[i];
            output.vals[r_out * output.ncols + c] = input.vals[i];
        }
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai, template<typename> typename Ap>
static void permute_matrix_cols(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, const Permutation<I,Ap> & perm)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != input.nvals) MY_ABORT("permute_matrix_cols: output matrix has wrong dimensions");

    // output[col] = input[perm[col]]

    struct colval{ I col; T val; colval(I c, T v){ col = c; val = v;} };

    std::copy(input.rowptrs, input.rowptrs + input.nrows + 1, output.rowptrs);

    for(I r = 0; r < output.nrows; r++)
    {
        I start = input.rowptrs[r];
        I end = input.rowptrs[r+1];
        I out_start = output.rowptrs[r];
        std::vector<colval> data;
        data.reserve(end - start);
        for(I i = start; i < end; i++) data.emplace_back(perm.backward[input.colidxs[i]], input.vals[i]);
        std::sort(data.begin(), data.end(), [&](const colval & l, const colval & r){ return l.col < r.col; });
        for(size_t j = 0; j < data.size(); j++) { output.colidxs[out_start + j] = data[j].col; output.vals[out_start + j] = data[j].val; }
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai, template<typename> typename Ap>
static void permute_matrix_rows(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, const Permutation<I,Ap> & perm, bool invperm = false)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != input.nvals) MY_ABORT("permute_matrix_rows: output matrix has wrong dimensions");

    // output[row] = input[perm[row]]

    const I * permvals;
    if(invperm) permvals = perm.backward;
    else        permvals = perm.forward;

    std::vector<I> out_row_nnz(output.nrows+1);
    for(I r_in = 0; r_in < input.nrows; r_in++)
    {
        I nnz = input.rowptrs[r_in+1] - input.rowptrs[r_in];
        out_row_nnz[permvals[r_in]] = nnz;
    }
    out_row_nnz[output.nrows] = 0;
    //std::exclusive_scan(out_row_nnz.begin(), out_row_nnz.end(), output.rowptrs, I{0});
    output.rowptrs[0] = 0;
    for(I r = 0; r < output.nrows; r++) output.rowptrs[r+1] = output.rowptrs[r] + out_row_nnz[r];

    for(I r_in = 0; r_in < input.nrows; r_in++)
    {
        I r_out = permvals[r_in];
        I instart = input.rowptrs[r_in];
        I inend = input.rowptrs[r_in+1];
        I outstart = output.rowptrs[r_out];
        std::copy(input.colidxs + instart, input.colidxs + inend, output.colidxs + outstart);
        std::copy(input.vals + instart, input.vals + inend, output.vals + outstart);
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai, template<typename> typename Ap>
static void permute_matrix(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, const Permutation<I,Ap> & perm)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != input.nvals) MY_ABORT("permute_matrix: output matrix has wrong dimensions");

    // probably does not even work, not tested

    struct colval{ I col; T val; colval(I c, T v){ col = c; val = v;} };

    I curr_idx = 0;
    for(I r_out = 0; r_out < output.nrows; r_out++)
    {
        output.rowptrs[r_out] = curr_idx;
        I r_in = perm.forward[r_out];
        I in_start = input.rowptrs[r_in];
        I in_end = input.rowptrs[r_in+1];
        std::vector<colval> colvals;
        colvals.reserve(in_end - in_start);
        for(I i = in_start; i < in_end; i++)
        {
            colvals.emplace_back(perm.backward(input.colidxs[i]), input.vals[i]);
        }
        std::sort(colvals.begin(), colvals.end(), [](const colval & l, const colval & r) { return l.col < r.col; });
        for(size_t j = 0; j < colvals.size(); j++)
        {
            output.colidxs[curr_idx] = colvals[j].col;
            output.vals[curr_idx] = colvals[j].val;
            curr_idx++;
        }
    }
    output.rowptrs[output.nrows] = curr_idx;
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai, template<typename> typename Ap>
static void permute_matrix_upper(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, const Permutation<I,Ap> & perm, putimers & tm)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != input.nvals) MY_ABORT("permute_matrix_upper: output matrix has wrong dimensions");
    struct colval{ I col; T val; colval(I c, T v){ col = c; val = v;} };

    tm.read.start();
    std::vector<std::vector<colval>> out_data(output.nrows);
    for(I r_in = 0; r_in < input.nrows; r_in++)
    {
        I r_out = perm.backward[r_in];
        I start = input.rowptrs[r_in];
        I end = input.rowptrs[r_in+1];
        for(I i = start; i < end; i++)
        {
            I c_in = input.colidxs[i];
            T val = input.vals[i];
            I c_out = perm.backward[c_in];
            if(r_out < c_out) out_data[r_out].emplace_back(c_out, val);
            else out_data[c_out].emplace_back(r_out, val);
        }
    }
    tm.read.stop();

    tm.sort.start();
    for(I r_out = 0; r_out < output.nrows; r_out++)
    {
        std::sort(out_data[r_out].begin(), out_data[r_out].end(), [&](const colval & l, const colval & r){ return l.col < r.col; });
    }
    tm.sort.stop();

    tm.write.start();
    I curr_idx = 0;
    for(I r_out = 0; r_out < output.nrows; r_out++)
    {
        output.rowptrs[r_out] = curr_idx;
        for(size_t j = 0; j < out_data[r_out].size(); j++)
        {
            output.colidxs[curr_idx] = out_data[r_out][j].col;
            output.vals[curr_idx] = out_data[r_out][j].val;
            curr_idx++;
        }
    }
    output.rowptrs[output.nrows] = curr_idx;
    tm.write.stop();
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai, template<typename> typename Ap>
static void permute_matrix_rows(MatrixDense<T,I,Ao> & output, const MatrixDense<T,I,Ai> & input, const Permutation<I,Ap> & perm, bool invperm = false)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.ld != input.ld) MY_ABORT("permute_matrix_rows: output matrix has wrong dimensions");
    // output[row] = input[perm[row]]

    const I * permvals;
    if(invperm) permvals = perm.backward;
    else        permvals = perm.forward;

    for(I r_in = 0; r_in < input.nrows; r_in++)
    {
        I r_out = permvals[r_in];
        std::copy(input.vals + r_in * input.ld, input.vals + r_in * input.ld + input.ncols, output.vals + r_out * output.ld);
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Aa, template<typename> typename Ab>
static void matrix_add(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Aa> & in1, const MatrixCSR<T,I,Ab> & in2)
{
    // add in2 to in1 and stores result in output, assumung all nonzero entries from in2 already exist in in1

    if(in1.nrows != in2.nrows || in1.ncols != in2.ncols) MY_ABORT("matrix_add: input matrices have non-matching dimensions");
    if(output.nrows != in1.nrows || output.ncols != in1.ncols || output.nvals != in1.nvals) MY_ABORT("matrix_add: output matrix has wrong dimensions");
    std::copy(in1.vals, in1.vals + in1.nvals, output.vals);
    std::copy(in1.colidxs, in1.colidxs + in1.nvals, output.colidxs);
    std::copy(in1.rowptrs, in1.rowptrs + in1.nrows + 1, output.rowptrs);

    for(I r = 0; r < in2.nrows; r++)
    {
        I i_start = in2.rowptrs[r];
        I i_end = in2.rowptrs[r+1];
        for(I i = i_start; i < i_end; i++)
        {
            I c = in2.colidxs[i];
            I * begin = output.colidxs + output.rowptrs[r];
            I * end = output.colidxs + output.rowptrs[r+1];
            I * res = std::find(begin, end, c);
            I idx = res - output.colidxs;
            if(res == end) MY_ABORT("matrix_add: in2 contains nz entry not present in in1");
            output.vals[idx] += in2.vals[i];
        }
    }
}



template<typename T, typename I, template<typename> typename A>
static void matrix_scale_inplace(MatrixDense<T,I,A> & matrix, T scalar)
{
    I nvals = matrix.size_alloc();
    for(I i = 0; i < nvals; i++)
    {
        matrix.vals[i] *= scalar;
    }
}



template<typename T, typename I, template<typename> typename Ai, template<typename> typename Ao>
static void matrix_scale_outofplace(MatrixDense<T,I,Ao> & output, const MatrixDense<T,I,Ai> & input, T scalar)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols) MY_ABORT("matrix_scale_outofplace: output matrix has wrong dimensions");

    if(output.ld == input.ld)
    {
        I nvals = input.size_alloc();
        for(I i = 0; i < nvals; i++)
        {
            output.vals[i] = scalar * input.vals[i];
        }
    }
    else
    {
        for(I r = 0; r < input.nrows; r++)
        {
            T * row_in = input.vals + r * input.ld;
            T * row_out = output.vals + r * output.ld;
            for(I c = 0; c < input.ncols; c++)
            {
                row_out[c] = scalar * row_in[c];
            }
        }
    }
}



template<typename I, template<typename> typename Ao, template<typename> typename Ai>
static void copy_permutation(Permutation<I,Ao> & output, const Permutation<I,Ai> & input)
{
    if(output.size != input.size) MY_ABORT("copy_permutation: output permutation has wrong dimension");
    std::copy(input.forward, input.forward + input.size, output.forward);
    std::copy(input.backward, input.backward + input.size, output.backward);
}



template<typename I, template<typename> typename Ao, template<typename> typename Ai>
static void inverse_permutation(Permutation<I,Ao> & output, const Permutation<I,Ai> & input)
{
    if(output.size != input.size) MY_ABORT("inverse_permutation: output permutation has wrong dimension");
    std::copy(input.forward, input.forward + input.size, output.backward);
    std::copy(input.backward, input.backward + input.size, output.forward);
}



template<typename I, template<typename> typename A>
static void inverse_permutation(Permutation<I,A> & perm, PermutationDirection todo)
{
    I *input, *output;
    if(todo == PermutationDirection::Forward)
    {
        input = perm.backward;
        output = perm.forward;
    }
    if(todo == PermutationDirection::Backward)
    {
        input = perm.forward;
        output = perm.backward;
    }

    for(I i = 0; i < perm.size; i++)
    {
        output[input[i]] = i;
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void matrix_transpose(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input)
{
    if(output.nrows != input.ncols || output.ncols != input.nrows || output.nvals != input.nvals) MY_ABORT("matrix_transpose: output matrix has wrong dimensions");
    struct colval{ I col; T val; colval(I c, T v){ col = c; val = v;} };

    std::vector<std::vector<colval>> out_rows(output.nrows);

    for(I r = 0; r < input.nrows; r++)
    {
        I start = input.rowptrs[r];
        I end = input.rowptrs[r+1];
        for(I i = start; i < end; i++)
        {
            I c = input.colidxs[i];
            T v = input.vals[i];
            out_rows[c].emplace_back(r, v);
        }
    }

    I curr_idx = 0;
    for(I row_out = 0; row_out < output.nrows; row_out++)
    {
        output.rowptrs[row_out] = curr_idx;
        std::vector<colval> & data = out_rows[row_out];
        for(size_t i = 0; i < data.size(); i++)
        {
            output.colidxs[curr_idx] = data[i].col;
            output.vals[curr_idx] = data[i].val;
            curr_idx++;
        }
    }
    output.rowptrs[output.nrows] = curr_idx;
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename A11, template<typename> typename A12, template<typename> typename A21, template<typename> typename A22>
static void concat_matrices(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,A11> & M11, const MatrixCSR<T,I,A12> & M12, const MatrixCSR<T,I,A21> & M21, const MatrixCSR<T,I,A22> & M22)
{
    if(M11.nrows != M12.nrows || M21.nrows != M22.nrows || M11.ncols != M21.ncols || M12.ncols != M22.ncols) MY_ABORT("concat_matrices: input matrices have wrong dimensions");
    if(output.nrows != M11.nrows + M21.nrows || output.ncols != M11.ncols + M12.ncols) MY_ABORT("concat_matrices: output matrix has wrong dimensions");

    I trows = M11.nrows;
    I brows = M21.nrows;
    I lcols = M11.ncols;
    I rcols = M12.ncols;
    I outnrows = trows + brows;
    I outncols = lcols + rcols;
    I outnvals = M11.nvals + M12.nvals + M21.nvals + M22.nvals;

    I curr_idx = 0;
    for(I r = 0; r < trows; r++)
    {
        output.rowptrs[r] = curr_idx;
        if(M11.nvals > 0)
        {
            I start = M11.rowptrs[r];
            I end = M11.rowptrs[r+1];
            for(I i = start; i < end; i++)
            {
                output.colidxs[curr_idx] = M11.colidxs[i];
                output.vals[curr_idx] = M11.vals[i];
                curr_idx++;
            }
        }
        if(M12.nvals > 0)
        {
            I start = M12.rowptrs[r];
            I end = M12.rowptrs[r+1];
            for(I i = start; i < end; i++)
            {
                output.colidxs[curr_idx] = M12.colidxs[i] + lcols;
                output.vals[curr_idx] = M12.vals[i];
                curr_idx++;
            }
        }
    }
    for(I r = 0; r < brows; r++)
    {
        I row_in = r;
        I row_out = r + trows;
        output.rowptrs[row_out] = curr_idx;
        if(M21.nvals > 0)
        {
            I start = M21.rowptrs[row_in];
            I end = M21.rowptrs[row_in+1];
            for(I i = start; i < end; i++)
            {
                output.colidxs[curr_idx] = M21.colidxs[i];
                output.vals[curr_idx] = M21.vals[i];
                curr_idx++;
            }
        }
        if(M22.nvals > 0)
        {
            I start = M22.rowptrs[row_in];
            I end = M22.rowptrs[row_in+1];
            for(I i = start; i < end; i++)
            {
                output.colidxs[curr_idx] = M22.colidxs[i] + lcols;
                output.vals[curr_idx] = M22.vals[i];
                curr_idx++;
            }
        }
    }
    output.rowptrs[output.nrows] = curr_idx;

    // {
    //     MatrixDense<T,I> d;
    //     sparse_to_dense(d, output);
    //     print_matrix(d, "concat output");
    // }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void submatrix(MatrixCSR<T,I,Ao> * output, const MatrixCSR<T,I,Ai> & input, I rowstart, I rowend, I colstart, I colend)
{
    if(output == nullptr) MY_ABORT("submatrix: output cannot be nullptr");

    if(rowstart < 0) rowstart = 0;
    if(rowend < 0) rowend = input.nrows;
    if(colstart < 0) colstart = 0;
    if(colend < 0) colend = input.ncols;

    I out_nrows = rowend - rowstart;
    I out_ncols = colend - colstart;
    std::vector<I> out_rowptrs(out_nrows + 1);
    std::vector<I> out_colidxs;
    std::vector<T> out_vals;

    I curr_nvals = 0;
    for(I r_out = 0; r_out < out_nrows; r_out++)
    {
        out_rowptrs[r_out] = curr_nvals;
        I r_in = r_out + rowstart;
        I start = input.rowptrs[r_in];
        I end = input.rowptrs[r_in+1];
        I i;
        for(i = start; i < end; i++) if(input.colidxs[i] >= colstart) break;
        for(; i < end; i++)
        {
            if(input.colidxs[i] >= colend) break;
            out_colidxs.push_back(input.colidxs[i] - colstart);
            out_vals.push_back(input.vals[i]);
            curr_nvals++;
        }
    }
    out_rowptrs[out_nrows] = curr_nvals;

    output->resize(out_nrows, out_ncols, curr_nvals, true);
    std::copy(out_rowptrs.begin(), out_rowptrs.end(), output->rowptrs);
    std::copy(out_colidxs.begin(), out_colidxs.end(), output->colidxs);
    std::copy(out_vals.begin(), out_vals.end(), output->vals);
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void submatrix(MatrixDense<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input, I rowstart, I rowend, I colstart, I colend)
{
    if(rowstart < 0) rowstart = 0;
    if(rowend < 0) rowend = input.nrows;
    if(colstart < 0) colstart = 0;
    if(colend < 0) colend = input.ncols;

    I out_nrows = rowend - rowstart;
    I out_ncols = colend - colstart;
    if(output.nrows != out_nrows || output.ncols != out_ncols) MY_ABORT("submatrix: output matrix has wrong dimensions");
    std::fill(output.vals, output.vals + output.size_alloc(), T{0});

    for(I r_out = 0; r_out < out_nrows; r_out++)
    {
        T * row_out = output.vals + r_out * output.ld;
        I r_in = r_out + rowstart;
        I start = input.rowptrs[r_in];
        I end = input.rowptrs[r_in+1];
        I i;
        for(i = start; i < end; i++) if(input.colidxs[i] >= colstart) break;
        for(; i < end; i++)
        {
            I c_in = input.colidxs[i];
            if(c_in >= colend) break;
            I c_out = c_in - colstart;
            T v = input.vals[i];
            row_out[c_out] = v;
        }
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void submatrix_select_rows(MatrixCSR<T,I,Ao> * output, const MatrixCSR<T,I,Ai> & input, const std::vector<I> & rows)
{
    // untested
    if(output == nullptr) MY_ABORT("submatrix: output cannot be nullptr");
    I out_nrows = rows.size();

    I output_nnz = 0;
    for(I r_out = 0; r_out < out_nrows; r_out++)
    {
        I r_in = rows[r_out];
        output_nnz += input.rowptrs[r_in+1] - input.rowptrs[r_in];
    }
    output->resize(out_nrows, input.ncols, output_nnz, true);

    I curr_idx = 0;
    for(I r_out = 0; r_out < out_nrows; r_out++)
    {
        output->rowptrs[r_out] = curr_idx;
        I r_in = rows[r_out];
        I start = input.rowptrs[r_in];
        I end = input.rowptrs[r_in+1];
        std::copy(input.colidxs + start, input.colidxs + end, output->colidxs + curr_idx);
        std::copy(input.vals + start, input.vals + end, output->vals + curr_idx);
        curr_idx += end - start;
    }
    output->rowptrs[output->nrows] = curr_idx;
}



template<typename T, typename I, template<typename> typename A>
static void upper_to_full(MatrixDense<T,I,A> & matrix)
{
    for(I r = 0; r < matrix.nrows; r++)
    {
        for(I c = 0; c < r; c++)
        {
            matrix.vals[r * matrix.ncols + c] = matrix.vals[c * matrix.ncols + r];
        }
    }
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void upper_to_full(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input)
{
    // assuming all diagonal nonzeros are present
    I diag_size = std::min(input.nrows, input.ncols);
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != 2 * input.nvals - diag_size) MY_ABORT("upper_to_full: output matrix has wrong dimensions");

    struct colval{ I col; T val; colval(I c, T v){ col = c; val = v;} };

    std::vector<std::vector<colval>> out_data(input.nrows);
    for(I r = 0; r < input.nrows; r++)
    {
        I start = input.rowptrs[r];
        I end = input.rowptrs[r+1];
        for(I i = start; i < end; i++)
        {
            I c = input.colidxs[i];
            T val = input.vals[i];
            out_data[r].emplace_back(c, val);
            if(r != c) out_data[c].emplace_back(r, val);
        }
    }

    // no need to sort, already sorted

    // I out_nvals = std::transform_reduce(out_data.begin(), out_data.end(), size_t{0}, std::plus<>(), [](const std::vector<colval> & v){return v.size();});
    I out_nvals = 0;
    for(size_t r = 0; r < out_data.size(); r++) out_nvals += out_data[r].size();
    output.resize(input.nrows, input.ncols, out_nvals, true);

    I curr_idx = 0;
    for(I r = 0; r < output.nrows; r++)
    {
        output.rowptrs[r] = curr_idx;
        for(size_t j = 0; j < out_data[r].size(); j++)
        {
            output.colidxs[curr_idx] = out_data[r][j].col;
            output.vals[curr_idx] = out_data[r][j].val;
            curr_idx++;
        }
    }
    output.rowptrs[output.nrows] = curr_idx;

    if(curr_idx != output.nvals) MY_ABORT("upper_to_full: input matrix probably does not contain all diagonal nonzeros");
}



template<typename T, typename I, template<typename> typename A>
static void print_matrix_csr(const MatrixCSR<T,I,A> & M, const char * name)
{
    static_assert(A<T>::is_data_host_accessible, "data have to be host accessible");
    static_assert(A<I>::is_data_host_accessible, "data have to be host accessible");
    printf("Matrix %s\n", name);
    printf("Size %lldx%lld, %lld nnz\n", (long long)M.nrows, (long long)M.ncols, (long long)M.nvals);
    for(I r = 0; r < M.nrows; r++)
    {
        printf("Row %4lld colidxs: ", (long long)r);
        for(I i = M.rowptrs[r]; i < M.rowptrs[r+1]; i++) printf("%6lld ", (long long)M.colidxs[i]);
        printf("\n");
        printf("Row %4lld vals   : ", (long long)r);
        for(I i = M.rowptrs[r]; i < M.rowptrs[r+1]; i++) printf("%+6.2lf ", (double)M.vals[i]);
        printf("\n");
    }
    fflush(stdout);
}



template<typename T, typename I, template<typename> typename A>
static void print_matrix_csr_arrays(const MatrixCSR<T,I,A> & M)
{
    // static_assert(A<T>::is_data_host_accessible, "data have to be host accessible");
    // static_assert(A<I>::is_data_host_accessible, "data have to be host accessible");

    printf("%lld %lld %lld", (long long)M.nrows, (long long)M.ncols, (long long)M.nvals);
    printf("\n");
    for(I r = 0; r <= M.nrows; r++) printf("%lld ", (long long)M.rowptrs[r]);
    printf("\n");
    for(I i = 0; i < M.nvals; i++) printf("%lld ", (long long)M.colidxs[i]);
    printf("\n");
    for(I i = 0; i < M.nvals; i++) printf("%.15e ", M.vals[i]);
    printf("\n");
    fflush(stdout);
}



template<typename T, typename I, template<typename> typename A>
static void print_matrix_ijv_arrays(const MatrixCSR<T,I,A> & M)
{
    // static_assert(A<T>::is_data_host_accessible, "data have to be host accessible");
    // static_assert(A<I>::is_data_host_accessible, "data have to be host accessible");

    printf("%lld %lld %lld", (long long)M.nrows, (long long)M.ncols, (long long)M.nvals);
    printf("\n");
    for(I r = 0; r < M.nrows; r++)
    {
        I start = M.rowptrs[r];
        I end = M.rowptrs[r+1];
        for(I i = start; i < end; i++) printf("%lld ", (long long)r);
    }
    printf("\n");
    for(I i = 0; i < M.nvals; i++) printf("%lld ", (long long)M.colidxs[i]);
    printf("\n");
    for(I i = 0; i < M.nvals; i++) printf("%.15e ", M.vals[i]);
    printf("\n");
    fflush(stdout);
}



template<typename T, typename I, template<typename> typename A>
static void print_matrix(const MatrixDense<T,I,A> & M, const char * name)
{
    static_assert(A<T>::is_data_host_accessible, "data have to be host accessible");
    printf("Matrix %s\n", name);
    printf("Size %lldx%lld, ld %lld\n", (long long)M.nrows, (long long)M.ncols, (long long)M.ld);
    for(I r = 0; r < M.nrows; r++)
    {
        for(I c = 0; c < M.ncols; c++)
        {
            if constexpr(std::is_floating_point_v<T>)
            {
                double v = (double)M.vals[r * M.ld + c];
                char str[100];
                sprintf(str, "%+10.3e", v);
                if(strstr(str, "nan") != nullptr) printf("  nan      ");
                else if(strstr(str, "inf") != nullptr) printf(" %cinf      ", v > 0 ? '+' : '-');
                else if(std::abs(v) >= 1e100) printf("  INF      ");
                else if(std::abs(v) <= 1e-100) printf("  0        ");
                else printf(" %+10.3e", v);
            }
            if constexpr(std::is_integral_v<T>) printf(" %+9lld", (long long)M.vals[r * M.ld + c]);
        }
        printf("\n");
    }
    fflush(stdout);
}



template<typename T, typename I, template<typename> typename A>
static void print_matrix(const MatrixCSR<T,I,A> & M, const char * name)
{
    static_assert(A<T>::is_data_host_accessible, "data have to be host accessible");
    MatrixDense<T,I,A> Mdense(M.nrows, M.ncols, I{-1}, true);
    sparse_to_dense(Mdense, M);
    print_matrix(Mdense, name);
}



template<typename I, template<typename> typename A>
static void print_permutation(const Permutation<I,A> & perm, const char * name)
{
    printf("Permutation %s, size %lld\n", name, (long long)perm.size);
    printf("Index:    "); for(I i = 0; i < perm.size; i++) printf(" %5lld", (long long)i); printf("\n");
    printf("Forward:  "); for(I i = 0; i < perm.size; i++) printf(" %5lld", (long long)perm.forward[i]); printf("\n");
    printf("Backward: "); for(I i = 0; i < perm.size; i++) printf(" %5lld", (long long)perm.backward[i]); printf("\n");
}



template<typename T, typename I, template<typename> typename A1, template<typename> typename A2>
static bool check_matrices_upper_equal(const MatrixDense<T,I,A1> & M1, const MatrixDense<T,I,A2> & M2, T tolerance, bool do_print_right, bool do_print_wrong)
{
    if(M1.nrows != M2.nrows || M1.ncols != M2.ncols)
    {
        if(do_print_wrong) printf("  Non-matching dimensions\n");
        return false;
    }

    T min_rel_err = std::numeric_limits<T>::max();
    T avg_rel_err = 0;
    T max_rel_err = 0;
    for(I r = 0; r < M1.nrows; r++)
    {
        T row_avg_rel_err = 0;
        for(I c = r; c < M1.ncols; c++)
        {
            T v1 = M1.vals[r * M1.ld + c];
            T v2 = M2.vals[r * M2.ld + c];
            T true_val = (v1 + v2) / 2;
            T abs_err = std::abs(v1 - v2);
            T rel_err = abs_err / std::abs(true_val);
            min_rel_err = std::min(min_rel_err, rel_err);
            row_avg_rel_err += rel_err;
            max_rel_err = std::max(max_rel_err, rel_err);
        }
        avg_rel_err += row_avg_rel_err;
    }
    avg_rel_err /= M1.nrows * M1.ncols;

    bool is_ok = (max_rel_err < tolerance);

    if(is_ok && do_print_right || !is_ok && do_print_wrong) printf("  Relative error: min %.3le, avg %.3le, max %.3le%s\n", (double)min_rel_err, (double)avg_rel_err, (double)max_rel_err, is_ok ? "" : " !!!");

    // {
    //     char str[100];
    //     sprintf(str, "%f", avg_rel_err);
    //     if(max_rel_err > 1e-10 || strstr(str, "nan") != nullptr)
    //     {
    //         printf("error too large, printing matrices\n");
    //         print_matrix(M1, "M1");
    //         print_matrix(M2, "M2");
    //     }
    // }

    return is_ok;
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void copy_matrix(MatrixCSR<T,I,Ao> & output, const MatrixCSR<T,I,Ai> & input)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != input.nvals) MY_ABORT("copy_matrix: output matrix has wrong dimensions");
    std::copy(input.rowptrs, input.rowptrs + input.nrows + 1, output.rowptrs);
    std::copy(input.colidxs, input.colidxs + input.nvals, output.colidxs);
    std::copy(input.vals, input.vals + input.nvals, output.vals);
}



template<typename To, typename Io, typename Ti, typename Ii, template<typename> typename Ao, template<typename> typename Ai>
static void cast_matrix(MatrixCSR<To,Io,Ao> & output, const MatrixCSR<Ti,Ii,Ai> & input)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols || output.nvals != input.nvals) MY_ABORT("cast_matrix: output matrix has wrong dimensions");
    std::transform(input.rowptrs, input.rowptrs + input.nrows + 1, output.rowptrs, [&](Ii r){ return static_cast<Io>(r); });
    std::transform(input.colidxs, input.colidxs + input.nvals, output.colidxs, [&](Ii c){ return static_cast<Io>(c); });
    std::transform(input.vals, input.vals + input.nvals, output.vals, [&](Ti v){ return static_cast<To>(v); });
}



template<typename T, typename I, template<typename> typename Ao, template<typename> typename Ai>
static void copy_matrix(MatrixDense<T,I,Ao> & output, const MatrixDense<T,I,Ai> & input)
{
    if(output.nrows != input.nrows || output.ncols != input.ncols) MY_ABORT("copy_matrix: output matrix has wrong dimensions");
    if(output.ld == input.ld) std::copy(input.vals, input.vals + input.size_alloc(), output.vals);
    else for(I r = 0; r < input.nrows; r++) std::copy(input.vals + r * input.ld, input.vals + r * input.ld + input.ncols, output.vals + r * output.ld);
}



static bool my_check_file_exists(const char * filepath)
{
    FILE * f = fopen(filepath, "r");
    bool ret;
    if(f == nullptr)
    {
        ret = false;
    }
    else
    {
        ret = true;
        fclose(f);
    }
    return ret;
}



template<typename T, typename I, template<typename> typename Ak, template<typename> typename Ab, template<typename> typename Af>
static void load_matrices(std::vector<MatrixCSR<T,I,Ak>> * out_Kregs, std::vector<MatrixCSR<T,I,Ab>> * out_Bs, std::vector<MatrixDense<T,I,Af>> * out_Fs_orig, const char * basematdir)
{
    if(out_Kregs == nullptr || out_Bs == nullptr || out_Fs_orig == nullptr) MY_ABORT("load_matrices: output vectors cannot be nullptr");

    MatrixCSR<T,I> K;
    MatrixCSR<T,I> RegMat;

    for(int i = 0; true; i++)
    {
        char pathB[100];
        char pathK[100];
        char pathRegMat[100];
        char pathF[100];

        snprintf(pathB, sizeof(pathB), "%s/B1%d.txt", basematdir, i);
        snprintf(pathK, sizeof(pathK), "%s/K%d.txt", basematdir, i);
        snprintf(pathRegMat, sizeof(pathRegMat), "%s/RegMat%d.txt", basematdir, i);
        snprintf(pathF, sizeof(pathF), "%s/F%d.txt", basematdir, i);

        if(!my_check_file_exists(pathK)) break;

        out_Bs->resize(out_Bs->size() + 1);
        out_Fs_orig->resize(out_Fs_orig->size() + 1);

        load_matrix_csr(&out_Bs->back(), pathB);
        load_matrix_csr(&K, pathK);
        load_matrix_csr(&RegMat, pathRegMat);
        load_matrix_dense(&out_Fs_orig->back(), pathF);

        out_Kregs->resize(out_Kregs->size() + 1);
        out_Kregs->back().resize(K.nrows, K.ncols, K.nvals, true);
        matrix_add(out_Kregs->back(), K, RegMat);
    }
}





template<typename T, typename I, template<typename> typename Ap, template<typename> typename Am>
static void my_sc_do_permutation(MatrixCSR<T,I,Ap> & Mperm, const MatrixCSR<T,I,Am> & M, I size_sc, permtimers & tm)
{
    if(Mperm.nrows != M.nrows || Mperm.ncols != M.ncols || Mperm.nvals != M.nvals) MY_ABORT("my_sc_do_permutation: output matrix has wrong dimensions");
    Permutation<I> perm(M.nrows, true);

    MatrixCSR<T,I> Mfull(M.nrows, M.ncols, 2 * M.nvals - M.nrows, true);
    tm.upper2full.start();
    upper_to_full(Mfull, M);
    tm.upper2full.stop();

    tm.vecinit.start();
    std::vector<bool> visited(Mfull.nrows, false);
    std::vector<bool> marked(Mfull.nrows, false);
    std::deque<I> que;
    std::vector<I> visited_order;
    visited_order.reserve(Mfull.nrows);
    tm.vecinit.stop();

    tm.traverse.start();
    I size_other = M.nrows - size_sc;
    for(I r = Mfull.nrows-1; r >= size_other; r--)
    {
        que.push_back(r);
        marked[r] = true;
    }

    while(!que.empty())
    {
        I r = que.front();
        que.pop_front();
        visited[r] = true;
        visited_order.push_back(r);
        I start = Mfull.rowptrs[r];
        I end = Mfull.rowptrs[r+1];
        for(I i = start; i < end; i++)
        {
            I c = Mfull.colidxs[i];
            if(!marked[c])
            {
                marked[c] = true;
                que.push_back(c);
            }
        }
    }
    tm.traverse.stop();

    tm.invperm.start();
    std::reverse_copy(visited_order.begin(), visited_order.end(), perm.forward);
    inverse_permutation(perm, PermutationDirection::Backward);
    tm.invperm.stop();

    tm.doperm.start();
    permute_matrix_upper(Mperm, M, perm, tm.pu);
    tm.doperm.stop();

    // {
    //     MatrixDense<T,I> m, mp;
    //     sparse_to_dense(m, M);
    //     sparse_to_dense(mp, Mperm);
    //     print_matrix(m, "M");
    //     print_matrix(mp, "Mperm");
    //     print_permutation(perm, "perm");
    // }
}



template<typename T, typename I, template<typename> typename Am>
static void my_sc_compute_boundaries(std::vector<I> * sc_boundary_idxs_, std::vector<I> * sc_sizes_, const MatrixCSR<T,I,Am> & Mperm, I sc_size, char scaling)
{
    if(sc_boundary_idxs_ == nullptr || sc_sizes_ == nullptr) MY_ABORT("my_sc_compute_boundaries: input vector cannot be nullptr");

    std::vector<I> & sc_boundary_idxs = *sc_boundary_idxs_;
    std::vector<I> & sc_sizes = *sc_sizes_;

    sc_boundary_idxs.reserve(Mperm.nrows / sc_size);
    sc_boundary_idxs.push_back(Mperm.nrows);
    sc_boundary_idxs.push_back(Mperm.nrows - sc_size);
    sc_sizes.reserve(Mperm.nrows / sc_size);
    sc_sizes.push_back(sc_size);
    while(sc_boundary_idxs.back() != 0)
    {
        I curr_sc_start = sc_boundary_idxs.back();
        I curr_sc_end = sc_boundary_idxs[sc_boundary_idxs.size() - 2];
        I curr_sc_size = curr_sc_end - curr_sc_start;
        I first_nz_row = 0;
        while(Mperm.colidxs[Mperm.rowptrs[first_nz_row+1]-1] < curr_sc_start) first_nz_row++;
        I min_inner_sc_size = curr_sc_start - first_nz_row;
        I inner_sc_size = -1;
        if(scaling == 'C') inner_sc_size = std::max(curr_sc_size, min_inner_sc_size);
        if(scaling == 'D') inner_sc_size = (curr_sc_start <= curr_sc_size) ? curr_sc_start : min_inner_sc_size;
        inner_sc_size = std::min(inner_sc_size, curr_sc_start);
        I inner_sc_end = curr_sc_start;
        I inner_sc_start = inner_sc_end - inner_sc_size;
        sc_boundary_idxs.push_back(inner_sc_start);
        sc_sizes.push_back(inner_sc_size);
    }
    std::reverse(sc_boundary_idxs.begin(), sc_boundary_idxs.end());
    std::reverse(sc_sizes.begin(), sc_sizes.end());

    // printf("SC boundaries:"); for(size_t i = 0; i < sc_boundary_idxs.size(); i++) printf(" %lld", (long long)sc_boundary_idxs[i]); printf("\n");
    // printf("SC sizes:"); for(size_t i = 0; i < sc_sizes.size(); i++) printf(" %lld", (long long)sc_sizes[i]); printf("\n");
}
