
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
static void applyB(FETI<T> &feti, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out, double * apply_start_time = nullptr)
{
    math::set(out, T{0});
    for (size_t d = 0; d < feti.K.size(); ++d) {
        for (int r = 0; r < feti.B1[d].nrows; ++r) {
            for (int c = feti.B1[d].rows[r]; c < feti.B1[d].rows[r + 1]; ++c) {
                out.vals[feti.D2C[d][r]] += feti.B1[d].vals[c] * in[d].vals[feti.B1[d].cols[c]];
            }
        }
    }
    if(apply_start_time != nullptr)
    {
        double stop = eslog::time();
        eslog::info("TMP DUAL OPERATOR APPLY TIME:  %12.6f ms\n", (stop - *apply_start_time) * 1000.0);
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
