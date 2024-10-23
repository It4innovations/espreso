
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
    indices.resize(feti.K.size());

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        Btx[d].resize(feti.K[d].nrows);
        KBtx[d].resize(feti.K[d].nrows);
        sc[d].shape = Matrix_Shape::UPPER;
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
    if (feti.updated.K || feti.updated.B) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            std::vector<int> permutation(feti.B1[d].ncols);
            for (int i = 0; i < feti.B1[d].nnz; ++i) {
                permutation[feti.B1[d].cols[i]] = 1;
            }
            indices[d].clear();
            for (size_t i = 0; i < permutation.size(); ++i) {
                if (permutation[i] == 1) {
                    indices[d].push_back(i);
                }
            }
            for (size_t i = 0, j = 0, k = feti.B1[d].ncols - indices[d].size(); i < permutation.size(); ++i) {
                permutation[i] = (permutation[i] == 1 ? k++ : j++);
            }
            Matrix_CSR<T> pK(feti.K[d]);
            math::permute(pK, permutation);

            sc[d].resize(indices[d].size(), indices[d].size());
            DirectSparseSolver<T> Ksolver;
            Ksolver.commit(pK);
            Ksolver.getSC(sc[d]);
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
        for (size_t i = 0; i < indices[d].size(); ++i) {
            Btx[d].vals[i] = Btx[d].vals[indices[d][i]];
        }
        math::blas::apply(KBtx[d], T{1}, sc[d], T{0}, Btx[d]);
        for (size_t i = indices[d].size(); i > 0; --i) {
            KBtx[d].vals[indices[d][i - 1]] = KBtx[d].vals[i - 1];
        }
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
