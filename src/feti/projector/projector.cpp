
#include "projector.h"
#include "orthogonal/tfetisymmetric.h"


namespace espreso {

template <typename T>
static Projector<T>* _set(FETI<T> *feti)
{
//	switch (feti->configuration.projector) {
		eslog::info(" = PROJECTOR                                                                      ORTHOGONAL = \n");
		return new OrthogonalTFETISymmetric<T>(feti);
//	}
}

template <> Projector<double>* Projector<double>::set(FETI<double> *feti) { return _set<double>(feti); }
template <> Projector<std::complex<double> >* Projector<std::complex<double> >::set(FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }

}
