
#include "regularization.h"
#include "esinfo/eslog.h"

namespace espreso {

template <typename T>
void Regularization<T>::set(step::Step &step)
{
	eslog::info(" = REGULARIZATION                                                                   ANALYTIC = \n");
	feti.R1.resize(feti.K.size());
	feti.R2.resize(feti.K.size());
	feti.RegMat.resize(feti.K.size());

	setAnalytic();
}

template <typename T>
void Regularization<T>::update(step::Step &step)
{
	updateAnalytic();
}

template <typename T>
void Regularization<T>::setAlgebraic()
{

}

template <typename T>
void Regularization<T>::updateAlgebraic()
{

}

template struct Regularization<double>;

}

