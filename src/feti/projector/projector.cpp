
#include "projector.h"
#include "orthogonal/tfetisymmetric.h"

namespace espreso {

template struct Projector<double>;
template struct Projector<std::complex<double> >;

template <typename T>
Projector<T>* Projector<T>::set(FETI<T> &feti, const step::Step &step)
{
//	switch (feti.configuration.projector) {
		eslog::info(" = PROJECTOR                                                                      ORTHOGONAL = \n");
		return new OrthogonalTFETISymmetric<T>(feti);
//	}
}

}
