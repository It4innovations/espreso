
#include "totalfeti.moore.penrose.h"
#include "feti/common/applyB.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "wrappers/mpi/communication.h"

namespace espreso {

template <typename T>
TotalFETIMoorePenrose<T>::TotalFETIMoorePenrose(FETI<T> &feti)
: DualOperator<T>(feti)
{

}

template <typename T>
TotalFETIMoorePenrose<T>::~TotalFETIMoorePenrose()
{

}

template <typename T>
void TotalFETIMoorePenrose<T>::info()
{
    if (this->infoPrinted && !feti.updated.K && !feti.updated.B) {
        return;
    }
    this->infoPrinted = true;

    eslog::info("      = ------------------------------------------------------------------------------- = \n");
    eslog::info("      = TOTAL FETI OPERATOR WITH MOORE-PENROSE INVERSE                                  = \n");
    eslog::info("      = ------------------------------------------------------------------------------- = \n");
}

/*
 * prepare buffers and call symbolic factorization that is independent on the Kplus values
 */
template <typename T>
void TotalFETIMoorePenrose<T>::set(const step::Step &step)
{
    Btx.resize(feti.K.size());
    KplusBtx.resize(feti.K.size());

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        Btx[di].resize(feti.K[di].nrows);
        KplusBtx[di].resize(feti.K[di].nrows);
        math::set(Btx[di], T{0});
    }

    eslog::checkpointln("FETI: SET TOTAL-FETI OPERATOR");
}

template <typename T>
void TotalFETIMoorePenrose<T>::update(const step::Step &step)
{
    if (feti.updated.B) {
        d.resize();
    }

    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{Kplus}\n");
        for (size_t di = 0; di < feti.K.size(); ++di) {
            math::store(feti.MoorePenroseInv[di], utils::filename(utils::debugDirectory(step) + "/feti/dualop", (std::string("Kplus") + std::to_string(di)).c_str()).c_str());
        }
    }

    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        math::blas::apply(KplusBtx[di], T{1}, feti.MoorePenroseInv[di], T{0}, feti.f[di]);
    }
    applyB(feti, KplusBtx, d);
    d.synchronize();
    math::add(d, T{-1}, feti.c);
    eslog::checkpointln("FETI: COMPUTE DUAL RHS [d]");
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/dualop/{d}\n");
        math::store(d, utils::filename(utils::debugDirectory(step) + "/feti/dualop", "d").c_str());
    }

    info();
}


template <typename T>
void TotalFETIMoorePenrose<T>::_apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
        math::blas::apply(KplusBtx[di], T{1}, feti.MoorePenroseInv[di], T{0}, Btx[di]);
    }
    applyB(feti, KplusBtx, y);
}

template <typename T>
void TotalFETIMoorePenrose<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    _apply(x, y);
    y.synchronize();
}

template <typename T>
void TotalFETIMoorePenrose<T>::apply(const Matrix_Dual<T> &x, Matrix_Dual<T> &y)
{
    Vector_Dual<T> _x, _y;
    _x.size = _y.size = x.ncols;
    for (int r = 0; r < x.nrows; ++r) {
        _x.vals = x.vals + x.ncols * r;
        _y.vals = y.vals + y.ncols * r;
        _apply(_x, _y);
    }
    y.synchronize();
}

template <typename T>
void TotalFETIMoorePenrose<T>::toPrimal(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y)
{
    #pragma omp parallel for
    for (size_t di = 0; di < feti.K.size(); ++di) {
        applyBt(feti, di, x, Btx[di]);
        math::copy(KplusBtx[di], feti.f[di]);
        math::add(KplusBtx[di], T{-1}, Btx[di]);
        math::blas::apply(y[di], T{1}, feti.MoorePenroseInv[di], T{0}, KplusBtx[di]);
    }
}

template class TotalFETIMoorePenrose<double>;

}
