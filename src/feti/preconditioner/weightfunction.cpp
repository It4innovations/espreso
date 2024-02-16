
#include "weightfunction.h"
#include "feti/common/applyB.h"

#include "esinfo/eslog.hpp"

namespace espreso {

template <typename T>
WeightFunction<T>::WeightFunction(FETI<T> &feti)
: Preconditioner<T>(feti)
{
    Btx.resize(feti.K.size());

    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        Btx[d].resize(feti.K[d].nrows);
    }

    eslog::checkpointln("FETI: SET WEIGHT FUNCTION PRECONDITIONER");
}

template <typename T>
void WeightFunction<T>::info()
{
    if (feti.configuration.exhaustive_info) {
        eslog::info(" = WEIGHT FUNCTION PRECONDITIONING PROPERTIES                                                = \n");
        eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    }
}

template <typename T>
void WeightFunction<T>::update(const step::Step &step)
{

}

template <typename T>
void WeightFunction<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
    #pragma omp parallel for
    for (size_t d = 0; d < feti.K.size(); ++d) {
        applyBt(feti, d, x, Btx[d]);
    }
    applyB(feti, Btx, y);
}

template struct WeightFunction<double>;
template struct WeightFunction<std::complex<double> >;

}
