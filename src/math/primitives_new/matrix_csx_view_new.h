
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_

#include <vector>

#include "math/primitives_new/matrix_base_new.h"
#include "math/primitives/matrix_csr.h"
#include "basis/utilities/utils.h"



namespace espreso {



template<typename T, typename I>
struct MatrixCsxView_new : public MatrixBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    I * ptrs = nullptr;
    I * idxs = nullptr;
    T * vals = nullptr;
    size_t nnz = 0;
    char order = '_';
    bool was_set = false;
public:
    using MatrixBase_new::nrows;
    using MatrixBase_new::ncols;
    using MatrixBase_new::prop;
public:
    MatrixCsxView_new() = default;
    MatrixCsxView_new(const MatrixCsxView_new &) = default;
    MatrixCsxView_new(MatrixCsxView_new &&) = default;
    MatrixCsxView_new & operator=(const MatrixCsxView_new &) = default;
    MatrixCsxView_new & operator=(MatrixCsxView_new &&) = default;
    virtual ~MatrixCsxView_new() = default;
public:
    void set_view(size_t nrows_, size_t ncols_, size_t nnz_, char order_, I * ptrs_, I * idxs_, T * vals_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized matrix view\n");
        nrows = nrows_;
        ncols = ncols_;
        nnz = nnz_;
        order = order_;
        ptrs = ptrs_;
        idxs = idxs_;
        vals = vals_;
        was_set = true;
    }
public:
    size_t get_size_primary() const
    {
        if(order == 'R') return nrows;
        if(order == 'C') return ncols;
        return 0;
    }
    size_t get_size_secdary() const
    {
        if(order == 'R') return ncols;
        if(order == 'C') return nrows;
        return 0;
    }

    MatrixCsxView_new<T,I> get_transposed_reordered_view() const
    {
        MatrixCsxView_new M = *this;
        M.transpose_reorder_inplace();
        return M;
    }
    void transpose_reorder_inplace()
    {
        std::swap(nrows, ncols);
        order = change_order(order);
        prop.uplo = change_uplo(prop.uplo);
    }

    static bool are_interchangable(const MatrixCsxView_new & A, const MatrixCsxView_new & B)
    {
        return (A.nrows == B.nrows) && (A.ncols == B.ncols) && (A.nnz == B.nnz) && (A.order == B.order) && (A.prop.uplo == B.prop.uplo) && (A.prop.diag == B.prop.diag);
    }

    template<typename A>
    static MatrixCsxView_new<T,I> from_old(const Matrix_CSR<T,I,A> & M_old)
    {
        MatrixCsxView_new<T,I> M_new;
        M_new.set_view(M_old.nrows, M_old.ncols, M_old.nnz, 'R', M_old.rows, M_old.cols, M_old.vals);
        if(M_old.shape == Matrix_Shape::LOWER) M_new.prop.uplo = 'L';
        if(M_old.shape == Matrix_Shape::UPPER) M_new.prop.uplo = 'U';
        if(M_old.shape == Matrix_Shape::FULL) M_new.prop.uplo = 'F';
        return M_new;
    }
    template<typename A>
    static Matrix_CSR<T,I,A> to_old(MatrixCsxView_new<T,I> & M_new)
    {
        Matrix_CSR<T,I,A> M_old;
        M_old.nrows = M_new.nrows;
        M_old.ncols = M_new.ncols;
        M_old.nnz = M_new.nnz;
        M_old.rows = M_new.ptrs;
        M_old.cols = M_new.idxs;
        M_old.vals = M_new.vals;
        if(M_new.prop.uplo == 'U') M_old.shape == Matrix_Shape::UPPER;
        else if(M_new.prop.uplo == 'L') M_old.shape == Matrix_Shape::LOWER;
        else M_old.shape == Matrix_Shape::FULL;
        return M_old;
    }

    void print(const char * name = "", char method = 'A')
    {
        if constexpr(utils::is_real<T>()) {
            eslog::info("CSX matrix %s, size %lldx%lld, nnz %lld, order '%c', uplo '%c', diag '%c'\n", name, nrows, ncols, nnz, order, prop.uplo, prop.diag);
            if(method == 'A') {
                eslog::info("ptrs: ");
                for(I ip = 0; ip <= get_size_primary(); ip++) eslog::info("%lld ", (long long)ptrs[ip]);
                eslog::info("\n");
                eslog::info("idxs: ");
                for(I i = 0; i < nnz; i++) eslog::info("%lld ", (long long)idxs[i]);
                eslog::info("\n");
                eslog::info("vals: ");
                for(I i = 0; i < nnz; i++) eslog::info("%+.3e ", (double)vals[i]);
                eslog::info("\n");
            }
            if(method == 'D') {
                struct rcv { I r; I c; T v; };
                std::vector<rcv> rcvs;
                for(I ip = 0; ip < get_size_primary(); ip++) {
                    I start = ptrs[ip];
                    I end = ptrs[ip+1];
                    for(I i = start; i < end; i++) {
                        I is = idxs[i];
                        T v = vals[i];
                        if(order == 'R') rcvs.push_back(rcv{ip,is,v});
                        if(order == 'C') rcvs.push_back(rcv{is,ip,v});
                    }
                }
                std::stable_sort(rcvs.begin(), rcvs.end(), [](const rcv & l, const rcv & r){ return l.c < r.c;});
                std::stable_sort(rcvs.begin(), rcvs.end(), [](const rcv & l, const rcv & r){ return l.r < r.r;});
                size_t curr_row = 0;
                size_t curr_col = 0;
                size_t curr_idx = 0;
                while(true) {
                    if(curr_idx < rcvs.size()) {
                        rcv & x = rcvs[curr_idx];
                        if(x.r == curr_row && x.c == curr_col) {
                            T v = x.v;
                            if(v == 0) eslog::info("    0       ");
                            else eslog::info(" %+11.3e", v);
                            curr_idx++;
                        }
                        else {
                            eslog::info("    .       ");
                        }
                    }
                    else {
                        eslog::info("    .       ");
                    }
                    curr_col++;
                    if(curr_col == ncols) {
                        curr_col = 0;
                        curr_row++;
                        eslog::info("\n");
                    }
                    if(curr_row == nrows) {
                        break;
                    }
                }
            }
            fflush(stdout);
        }
        if constexpr(utils::is_complex<T>()) {
            eslog::error("matrix print not supported for complex matrices\n");
        }
    }
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_ */
