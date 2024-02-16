
#include "dirichlet.h"
#include "feti/common/applyB.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/utilities/sysutils.h"

#include <algorithm>

namespace espreso {

template <typename T>
Dirichlet<T>::Dirichlet(FETI<T> &feti)
: Preconditioner<T>(feti)
{
    if (feti.configuration.ordering == FETIConfiguration::ORDERING::NATURAL) {
        eslog::error("natural ordering is not compatible with the Dirichlet preconditioner.\n");
    }
    Btx.resize(feti.K.size());
    KBtx.resize(feti.K.size());
    sc.resize(feti.K.size());
    Ksolver.resize(feti.K.size());

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        Ksolver[d].commit(feti.K[d]);
        Btx[d].resize(feti.K[d].nrows);
        KBtx[d].resize(feti.K[d].nrows);
        sc[d].shape = Matrix_Shape::UPPER;

        esint sc_size = feti.B1[d].ncols - *std::min_element(feti.B1[d].cols, feti.B1[d].cols + feti.B1[d].nnz);
        sc[d].resize(sc_size, sc_size);
    }

    eslog::checkpointln("FETI: SET DIRICHLET PRECONDITIONER");
}

template <typename T>
Dirichlet<T>::~Dirichlet()
{

}

template <typename T>
void Dirichlet<T>::info()
{
    if (feti.configuration.exhaustive_info) {
        eslog::info(" = DIRICHLET PRECONDITIONER PROPERTIES                                                       = \n");
        eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    }
}

template <typename T>
void Dirichlet<T>::update(const step::Step &step)
{
    if (feti.updated.K) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            Ksolver[d].getSC(sc[d]);
        }
    }
    _print(step);
}

template <typename T>
void Dirichlet<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        applyBt(feti, d, x, Btx[d]);
        Vector_Dense<T> _y, _x;
        _y.vals = KBtx[d].vals + KBtx[d].size - sc[d].nrows;
        _x.vals = Btx[d].vals + Btx[d].size - sc[d].nrows;
        math::blas::apply(_y, T{1}, sc[d], T{0}, _x);
    }
    applyB(feti, KBtx, y);
}

template <typename T>
void Dirichlet<T>::_print(const step::Step &step)
{
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: feti/preconditioner/{Dirichlet}\n");
        for (size_t d = 0; d < feti.K.size(); ++d) {
            math::store(sc[d], utils::filename(utils::debugDirectory(step) + "/feti/precondition", (std::string("dirichlet") + std::to_string(d)).c_str()).c_str());
        }
    }
}

template struct Dirichlet<double>;
template struct Dirichlet<std::complex<double> >;

}
