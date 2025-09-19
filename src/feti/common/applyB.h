
#ifndef SRC_FETI_COMMON_APPLYB_H_
#define SRC_FETI_COMMON_APPLYB_H_

#include "vector_dual.h"
#include "feti/feti.h"
#include "math/math.h"

namespace espreso {

template <typename T>
static void applyBt(FETI<T> &feti, size_t d, const Vector_Dual<T> &in, Vector_Dense<T> &out, T alpha = T{1})
{
    math::set(out, T{0});

    for (int r = 0; r < feti.B1[d].nrows; ++r) {
        for (int c = feti.B1[d].rows[r]; c < feti.B1[d].rows[r + 1]; ++c) {
            out.vals[feti.B1[d].cols[c]] += alpha * feti.B1[d].vals[c] * in.vals[feti.D2C[d][r]];
        }
    }
}

template <typename T>
static void extractDomain(FETI<T> &feti, size_t d, const Vector_Dual<T> &in, Vector_Dense<T> &out)
{
    math::set(out, T{0});
    for (int r = 0; r < feti.B1[d].nrows; ++r) {
        out.vals[r] = in.vals[feti.D2C[d][r]];
    }
}

// TODO: threaded implementation + more efficient 'beta' scale
template <typename T>
static void applyB(FETI<T> &feti, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out)
{
    math::set(out, T{0});
    for (size_t d = 0; d < feti.K.size(); ++d) {
        for (int r = 0; r < feti.B1[d].nrows; ++r) {
            for (int c = feti.B1[d].rows[r]; c < feti.B1[d].rows[r + 1]; ++c) {
                out.vals[feti.D2C[d][r]] += feti.B1[d].vals[c] * in[d].vals[feti.B1[d].cols[c]];
            }
        }
    }
}

template <typename T>
static void applyB(FETI<T> &feti, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out, const std::vector<int> &filter)
{
    math::set(out, T{0});
    for (size_t di = 0; di < filter.size(); ++di) {
        for (int r = 0; r < feti.B1[filter[di]].nrows; ++r) {
            for (int c = feti.B1[filter[di]].rows[r]; c < feti.B1[filter[di]].rows[r + 1]; ++c) {
                out.vals[feti.D2C[filter[di]][r]] += feti.B1[filter[di]].vals[c] * in[filter[di]].vals[feti.B1[filter[di]].cols[c]];
            }
        }
    }
}

template <typename T>
static void applyB_threaded(FETI<T> &feti, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out)
{
    math::set(out, T{0});
    #pragma omp parallel for schedule(static,1)
    for (size_t d = 0; d < feti.K.size(); ++d) {
        for (int r = 0; r < feti.B1[d].nrows; ++r) {
            for (int c = feti.B1[d].rows[r]; c < feti.B1[d].rows[r + 1]; ++c) {
                T & dst = out.vals[feti.D2C[d][r]];
                T val = feti.B1[d].vals[c] * in[d].vals[feti.B1[d].cols[c]];
                #pragma omp atomic
                dst += val;
            }
        }
    }
}

template <typename T>
static void insertDomains(FETI<T> &feti, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out)
{
    math::set(out, T{0});
    for (size_t d = 0; d < feti.K.size(); ++d) {
        for (int r = 0; r < feti.B1[d].nrows; ++r) {
            out.vals[feti.D2C[d][r]] += in[d].vals[r];
        }
    }
}

}

#endif /* SRC_FETI_COMMON_APPLYB_H_ */
