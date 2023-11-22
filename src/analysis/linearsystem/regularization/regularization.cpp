
#include "regularization.h"

namespace espreso {

template <typename T>
void Regularization<T>::set(step::Step &step)
{
	eslog::info(" = REGULARIZATION                                                                   ANALYTIC = \n");
	feti.regularization.R1.domains.resize(feti.K.domains.size());
	feti.regularization.R2.domains.resize(feti.K.domains.size());
	feti.regularization.RegMat.domains.resize(feti.K.domains.size());

	feti.regularization.RegMat.type = feti.K.type;
	feti.regularization.RegMat.shape = feti.K.shape;

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

