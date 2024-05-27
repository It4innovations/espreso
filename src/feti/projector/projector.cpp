
#include "projector.h"
#include "hfeti.orthogonal.symmetric.h"
#include "tfeti.orthogonal.symmetric.h"
#include "basis/utilities/sysutils.h"
#include "basis/utilities/utils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "math/math.h"

namespace espreso {

template <typename T> int Projector<T>::Kernel::roffset = 0;
template <typename T> int Projector<T>::Kernel::rsize   = 0;
template <typename T> int Projector<T>::Kernel::total   = 0;

template <typename T>
Projector<T>* Projector<T>::set(FETI<T> &feti, const step::Step &step)
{
    switch (feti.configuration.method) {
    case FETIConfiguration::METHOD::TOTAL_FETI: {
        switch (feti.configuration.projector) {
        case FETIConfiguration::PROJECTOR::ORTHOGONAL:
            eslog::info(" = PROJECTOR                                                             EXPLICIT ORTHOGONAL = \n");
            return new TFETIOrthogonalSymmetric<T>(feti);
        default: return nullptr;
        }
    } break;
    case FETIConfiguration::METHOD::HYBRID_FETI: {
        switch (feti.configuration.projector) {
        case FETIConfiguration::PROJECTOR::ORTHOGONAL:
            if (feti.configuration.projector_opt & FETIConfiguration::PROJECTOR_OPT::FULL) {
                eslog::info(" = PROJECTOR                                                        FULL EXPLICIT ORTHOGONAL = \n");
                return new TFETIOrthogonalSymmetric<T>(feti);
            } else {
                eslog::info(" = PROJECTOR                                                      HYBRID EXPLICIT ORTHOGONAL = \n");
                return new HFETIOrthogonalSymmetric<T>(feti);
            }
        default: return nullptr;
        }
    } break;
    }
    return nullptr;
}

template<typename T>
Projector<T>::Projector(FETI<T> &feti)
: feti(feti), GGTtype(GGT_TYPE::GGT)
{
    if (feti.configuration.projector_opt & FETIConfiguration::PROJECTOR_OPT::WITH_FACTORS) {
        GGTtype = GGT_TYPE::LU;
    }
    if (feti.configuration.iterative_solver == FETIConfiguration::ITERATIVE_SOLVER::SMALBE) {
        GGTtype = GGT_TYPE::LU;
    }
}

template<typename T>
void Projector<T>::info()
{
    eslog::info(" = PROJECTOR PROPERTIES                                                                      = \n");
    eslog::info(" =   ROWS                                                                          %9d = \n", std::max(invGGt.ncols, invL.ncols));
    eslog::info(" =   NNZ                                                                           %9d = \n", GGt.nnz);
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
}

template<typename T>
void Projector<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    x.copyToWithoutHalo(y);
    _apply_G(x, Gx);
    _apply_invGGt(Gx, iGGtGx);
    _apply_Gt(iGGtGx, T{-1}, y);
}

template<typename T>
void Projector<T>::apply_GtintGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
    math::set(y, T{0});
    _apply_invGGt(x, iGGtGx);
    _apply_Gt(iGGtGx, T{1}, y);
}

template<typename T>
void Projector<T>::apply_R(const Vector_Kernel<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if (GGTtype != GGT_TYPE::NONE) { // avoid this condition
        _apply_R(x, y);
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < y.size(); ++d) {
            math::set(y[d], T{0});
        }
    }
}

template<typename T>
void Projector<T>::apply_RinvGGtG(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    if (GGTtype != GGT_TYPE::NONE) { // avoid this condition
        _apply_G(x, Gx);
        _apply_invGGt(Gx, iGGtGx);
        _apply_R(iGGtGx, y);
    } else {
        #pragma omp parallel for
        for (size_t d = 0; d < y.size(); ++d) {
            math::set(y[d], T{0});
        }
    }
}

template<typename T>
void Projector<T>::apply_invU(const Vector_Kernel<T> &x, Vector_Kernel<T> &y)
{
    math::set(y, T{0});
    Vector_Dense<T> local;
    local.size = Vector_Kernel<T>::localSize;
    local.vals = y.vals + Vector_Kernel<T>::offset;
    _apply_invU(x, local);
    y.synchronize();
}

template<typename T>
void Projector<T>::apply_invL(const Vector_Kernel<T> &x, Vector_Kernel<T> &y)
{
    math::set(y, T{0});
    Vector_Dense<T> local;
    local.size = Vector_Kernel<T>::localSize;
    local.vals = y.vals + Vector_Kernel<T>::offset;
    _apply_invL(x, local);
    y.synchronize();
}

template<typename T>
void Projector<T>::apply_GtinvU(const Vector_Kernel<T> &x, Vector_Dual<T> &y)
{
    math::set(y, T{0});
    _apply_invU(x, iGGtGx);
    _apply_Gt(iGGtGx, T{1}, y);
}

template<typename T>
void Projector<T>::apply_invLG(const Vector_Dual<T> &x, Vector_Kernel<T> &y)
{
    _apply_G(x, Gx);
    apply_invL(Gx, y);
}

template<typename T>
void Projector<T>::_apply_G(const Vector_Dual<T> &in, Vector_Kernel<T> &out)
{
    #pragma omp parallel for
    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t r = Vector_Kernel<T>::distribution[t]; r < Vector_Kernel<T>::distribution[t + 1]; ++r) {
            out.vals[r + Vector_Kernel<T>::offset] = T{0};
            for (int c = G.rows[r]; c < G.rows[r + 1]; ++c) {
                out.vals[r + Vector_Kernel<T>::offset] += G.vals[c] * in.vals[G.cols[c]];
            }
        }
    }
    out.synchronize();
}

template<typename T>
void Projector<T>::_apply_invGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
    switch (GGTtype) {
    case GGT_TYPE::GGT:
        #pragma omp parallel for
        for (int t = 0; t < info::env::threads; ++t) {
            Matrix_Dense<T> a;
            Vector_Dense<T> y;
            a.ncols = invGGt.ncols;
            a.nrows = y.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];

            a.vals = invGGt.vals + invGGt.ncols * Vector_Kernel<T>::distribution[t];
            y.vals = out.vals + Vector_Kernel<T>::distribution[t];

            math::blas::apply(y, T{1}, a, T{0}, in);
        }
        break;
    case GGT_TYPE::LU: {
        Vector_Kernel<T> mid;
        #pragma omp parallel for
        for (int t = 0; t < info::env::threads; ++t) {
            Matrix_Dense<T> a;
            Vector_Dense<T> y;
            a.ncols = invL.ncols;
            a.nrows = y.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];
            y.vals = mid.vals + Vector_Kernel<T>::distribution[t] + Vector_Kernel<T>::offset;

            a.vals = invL.vals + invL.ncols * Vector_Kernel<T>::distribution[t];
            math::blas::apply(y, T{1}, a, T{0}, in);
        }
        mid.synchronize();
        #pragma omp parallel for
        for (int t = 0; t < info::env::threads; ++t) {
            Matrix_Dense<T> a;
            Vector_Dense<T> x;
            a.ncols = invU.ncols;
            a.nrows = x.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];
            x.vals = out.vals + Vector_Kernel<T>::distribution[t];

            a.vals = invU.vals + invU.ncols * Vector_Kernel<T>::distribution[t];
            math::blas::apply(x, T{1}, a, T{0}, mid);
        }
        } break;
    case GGT_TYPE::NONE:
        break;
    }
}

template<typename T>
void Projector<T>::_apply_invL(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
    switch (GGTtype) {
    case GGT_TYPE::LU:
        #pragma omp parallel for
        for (int t = 0; t < info::env::threads; ++t) {
            Matrix_Dense<T> a;
            Vector_Dense<T> x;
            a.ncols = invL.ncols;
            a.nrows = x.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];
            x.vals = out.vals + Vector_Kernel<T>::distribution[t];

            a.vals = invL.vals + invL.ncols * Vector_Kernel<T>::distribution[t];
            math::blas::apply(x, T{1}, a, T{0}, in);
        }
        break;
    case GGT_TYPE::GGT:
        eslog::error("cannot apply inv(L).\n");
        break;
    case GGT_TYPE::NONE:
        break;
    }
}

template<typename T>
void Projector<T>::_apply_invU(const Vector_Kernel<T> &in, Vector_Dense<T> &out)
{
    switch (GGTtype) {
    case GGT_TYPE::LU:
        #pragma omp parallel for
        for (int t = 0; t < info::env::threads; ++t) {
            Matrix_Dense<T> a;
            Vector_Dense<T> x;
            a.ncols = invU.ncols;
            a.nrows = x.size = Vector_Kernel<T>::distribution[t + 1] - Vector_Kernel<T>::distribution[t];
            x.vals = out.vals + Vector_Kernel<T>::distribution[t];

            a.vals = invU.vals + invU.ncols * Vector_Kernel<T>::distribution[t];
            math::blas::apply(x, T{1}, a, T{0}, in);
        }
        break;
    case GGT_TYPE::GGT:
        eslog::error("cannot apply inv(U).\n");
        break;
    case GGT_TYPE::NONE:
        break;
    }
}

template<typename T>
void Projector<T>::_apply_Gt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out)
{
    for (int r = 0; r < G.nrows; ++r) {
        for (int c = G.rows[r]; c < G.rows[r + 1]; ++c) {
            out.vals[G.cols[c]] += alpha * G.vals[c] * in.vals[r];
        }
    }
    out.synchronize();
}

template<typename T>
void Projector<T>::_apply_R(const Vector_Dense<T> &in, std::vector<Vector_Dense<T> > &out)
{
    #pragma omp parallel for
    for (size_t d = 0; d < out.size(); ++d) {
        Vector_Dense<T> y;
        y.size = feti.R1[d].nrows;
        y.vals = in.vals + Kernel::roffset + kernel[d].offset;

        math::blas::applyT(out[d], T{1}, feti.R1[d], T{0}, y);
    }
}

template<typename T>
void Projector<T>::_print(const step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/projector/{G, e, GGt, invGGt}\n");
        math::store(G, utils::filename(utils::debugDirectory(step) + "/feti/projector", "G").c_str());
        math::store(Gt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "Gt").c_str());
        math::store(e, utils::filename(utils::debugDirectory(step) + "/feti/projector", "e").c_str());
        math::store(GGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "GGt").c_str());
        switch (Projector<T>::GGTtype) {
        case Projector<T>::GGT_TYPE::GGT:
            math::store(invGGt, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invGGt").c_str());
            break;
        case Projector<T>::GGT_TYPE::LU:
            math::store(invU, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invU").c_str());
            math::store(invL, utils::filename(utils::debugDirectory(step) + "/feti/projector", "invL").c_str());
            break;
        default: break;
        }
    }
}

template struct Projector<double>;
template struct Projector<std::complex<double> >;

}
