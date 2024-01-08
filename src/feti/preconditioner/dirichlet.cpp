
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
	Btx.resize(feti.K.domains.size());
	KBtx.resize(feti.K.domains.size());
	sc.resize(feti.K.domains.size());
	Ksolver.resize(feti.K.domains.size());

	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		Ksolver[d].commit(feti.K.domains[d]);
		Btx[d].resize(feti.K.domains[d].nrows);
		KBtx[d].resize(feti.K.domains[d].nrows);
		sc[d].shape = Matrix_Shape::UPPER;

		const typename FETI<T>::EqualityConstraints::Domain &L = feti.equalityConstraints.domain[d];
		esint sc_size = L.B1.ncols - *std::min_element(L.B1.cols, L.B1.cols + L.B1.nnz);
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
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
		Ksolver[d].getSC(sc[d]);
	}
	_print(step);
}

template <typename T>
void Dirichlet<T>::apply(const Vector_Dual<T> &x, Vector_Dual<T> &y)
{
	#pragma omp parallel for
	for (size_t d = 0; d < feti.K.domains.size(); ++d) {
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
		for (size_t d = 0; d < feti.K.domains.size(); ++d) {
			math::store(sc[d], utils::filename(utils::debugDirectory(step) + "/feti/precondition", (std::string("Dirichlet") + std::to_string(d)).c_str()).c_str());
		}
	}
}

template struct Dirichlet<double>;
template struct Dirichlet<std::complex<double> >;

}
