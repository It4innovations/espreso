
#ifndef SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_MATRIX_CSX_VIEW_NEW_H_

#include <vector>

#include "math/primitives_new/matrix_base_new.h"
#include "math/primitives_new/allocator_new_base.h"
#include "math/primitives/matrix_csr.h"
#include "basis/utilities/utils.h"



namespace espreso {



template<typename T, typename I>
struct MatrixCsxView_new : public MatrixBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    Allocator_new * ator = nullptr;
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
    void set_view(size_t nrows_, size_t ncols_, size_t nnz_, char order_, I * ptrs_, I * idxs_, T * vals_, Allocator_new * ator_)
    {
        if(was_set) eslog::error("can only set yet-uninitialized matrix view\n");
        nrows = nrows_;
        ncols = ncols_;
        nnz = nnz_;
        order = order_;
        ptrs = ptrs_;
        idxs = idxs_;
        vals = vals_;
        ator = ator_;
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

    MatrixCsxView_new<T,I> get_transposed_reordered_view(bool do_conj = false) const
    {
        MatrixCsxView_new M = *this;
        M.transpose_reorder_inplace(do_conj);
        return M;
    }
    void transpose_reorder_inplace(bool do_conj = false)
    {
        std::swap(nrows, ncols);
        order = change_order(order);
        conj = (conj != do_conj);
        prop.uplo = change_uplo(prop.uplo);
    }
public:
    template<typename A>
    static MatrixCsxView_new<T,I> from_old(const Matrix_CSR<T,I,A> & M_old)
    {
        MatrixCsxView_new<T,I> M_new;
        M_new.set_view(M_old.nrows, M_old.ncols, M_old.nnz, 'R', M_old.rows, M_old.cols, M_old.vals, AllocatorDummy_new::get_singleton(A::is_data_host_accessible, A::is_data_device_accessible));
        M_new.prop.uplo = get_new_matrix_uplo(M_old.shape);
        M_new.prop.diag = '_';
        M_new.prop.symm = get_new_matrix_symmetry(M_old.type);
        M_new.prop.dfnt = get_new_matrix_definitness(M_old.type);
        return M_new;
    }
    template<typename A>
    static Matrix_CSR<T,I,A> to_old(const MatrixCsxView_new<T,I> & M_new)
    {
        if(M_new.order != 'R') eslog::error("can only convert to old row-major matrices\n");
        if(A::is_data_host_accessible != M_new.ator->is_data_accessible_cpu()) eslog::error("allocator access mismatch on cpu\n");
        if(A::is_data_device_accessible != M_new.ator->is_data_accessible_gpu()) eslog::error("allocator access mismatch on gpu\n");
        Matrix_CSR<T,I,A> M_old;
        M_old._allocated.nrows = M_new.nrows;
        M_old._allocated.ncols = M_new.ncols;
        M_old._allocated.nnz = M_new.nnz;
        M_old.nrows = M_new.nrows;
        M_old.ncols = M_new.ncols;
        M_old.nnz = M_new.nnz;
        M_old.rows = M_new.ptrs;
        M_old.cols = M_new.idxs;
        M_old.vals = M_new.vals;
        M_old.shape = get_old_matrix_shape(M_new);
        M_old.type = get_old_matrix_type<T>(M_new);
        return M_old;
    }
public:
    void print(const char * name = "", char method = 'A')
    {
        if(!ator->is_data_accessible_cpu()) eslog::error("print is supported only for cpu-accessible matrices\n");
        if constexpr(utils::is_real<T>()) {
            eslog::info("CSX matrix %s, size %lldx%lld, nnz %lld, order '%c', uplo '%c', diag '%c'\n", name, nrows, ncols, nnz, order, prop.uplo, prop.diag);
            if(method == 'A') {
                eslog::info("ptrs: ");
                for(size_t ip = 0; ip <= get_size_primary(); ip++) eslog::info("%lld ", (long long)ptrs[ip]);
                eslog::info("\n");
                eslog::info("idxs: ");
                for(size_t i = 0; i < nnz; i++) eslog::info("%lld ", (long long)idxs[i]);
                eslog::info("\n");
                eslog::info("vals: ");
                for(size_t i = 0; i < nnz; i++) eslog::info("%+.3e ", (double)vals[i]);
                eslog::info("\n");
            }
            if(method == 'D' || method == 'P') {
                struct rcv { I r; I c; T v; };
                std::vector<rcv> rcvs;
                for(size_t ip = 0; ip < get_size_primary(); ip++) {
                    I start = ptrs[ip];
                    I end = ptrs[ip+1];
                    for(I i = start; i < end; i++) {
                        I is = idxs[i];
                        T v = vals[i];
                        if(order == 'R') rcvs.push_back(rcv{(I)ip,is,v});
                        if(order == 'C') rcvs.push_back(rcv{is,(I)ip,v});
                    }
                }
                std::stable_sort(rcvs.begin(), rcvs.end(), [](const rcv & l, const rcv & r){ return l.c < r.c;});
                std::stable_sort(rcvs.begin(), rcvs.end(), [](const rcv & l, const rcv & r){ return l.r < r.r;});
                I curr_row = 0;
                I curr_col = 0;
                size_t curr_idx = 0;
                while(true) {
                    if(curr_idx < rcvs.size()) {
                        rcv & x = rcvs[curr_idx];
                        if(x.r == curr_row && x.c == curr_col) {
                            T v = x.v;
                            if(method == 'D') {
                                if(v == 0) eslog::info("    0       ");
                                else eslog::info(" %+11.3e", v);
                            }
                            else {
                                eslog::info("X");
                            }
                            curr_idx++;
                        }
                        else {
                            if(method == 'D') eslog::info("    .       ");
                            if(method == 'P') eslog::info(".");
                        }
                    }
                    else {
                        if(method == 'D') eslog::info("    .       ");
                            if(method == 'P') eslog::info(".");
                    }
                    curr_col++;
                    if(curr_col == (I)ncols) {
                        curr_col = 0;
                        curr_row++;
                        eslog::info("\n");
                    }
                    if(curr_row == (I)nrows) {
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
