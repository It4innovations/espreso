
#include "projector.h"
#include "orthogonal/tfetisymmetric.h"


namespace espreso {

template <typename T>
static Projector<T>* _set(AX_FETI<T> *feti)
{
//	switch (feti->configuration.projector) {
		eslog::info(" = PROJECTOR                                                                      ORTHOGONAL = \n");
		return new OrthogonalTFETISymmetric<T>(feti);
//	}
}

template <> Projector<double>* Projector<double>::set(AX_FETI<double> *feti) { return _set<double>(feti); }
template <> Projector<std::complex<double> >* Projector<std::complex<double> >::set(AX_FETI<std::complex<double> > *feti) { return _set<std::complex<double> >(feti); }

}
