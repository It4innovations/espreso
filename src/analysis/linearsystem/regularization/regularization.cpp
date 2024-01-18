
#include "regularization.h"

namespace espreso {

template <typename T>
void Regularization<T>::set(step::Step &step)
{
	eslog::info(" = REGULARIZATION                                                                   ANALYTIC = \n");
	feti.regularization.R1.resize(feti.K.domains.size());
	feti.regularization.R2.resize(feti.K.domains.size());
	feti.regularization.RegMat.resize(feti.K.domains.size());

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

