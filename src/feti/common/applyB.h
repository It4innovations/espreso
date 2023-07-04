
#ifndef SRC_FETI_COMMON_APPLYB_H_
#define SRC_FETI_COMMON_APPLYB_H_

#include "feti/feti.h"
#include "math/feti/vector_dual.h"

namespace espreso {

template <typename T>
static void applyBt(FETI<T> &feti, size_t d, const Vector_Dual<T> &in, Vector_Dense<T> &out)
{
	const typename FETI<T>::EqualityConstraints::Domain &L = feti.equalityConstraints.domain[d];

	math::set(out, T{0});

	for (esint r = 0; r < L.B1.nrows; ++r) {
		for (esint c = L.B1.rows[r]; c < L.B1.rows[r + 1]; ++c) {
			out.vals[L.B1.cols[c]] += L.B1.vals[c] * in.vals[L.D2C[r]];
		}
	}
}

template <typename T>
static void extractDomain(FETI<T> &feti, size_t d, const Vector_Dual<T> &in, Vector_Dense<T> &out)
{
	const typename FETI<T>::EqualityConstraints::Domain &L = feti.equalityConstraints.domain[d];

	math::set(out, T{0});
	for (esint r = 0; r < L.B1.nrows; ++r) {
		out.vals[r] = in.vals[L.D2C[r]];
	}
}

// TODO: threaded implementation + more efficient 'beta' scale
template <typename T>
static void applyB(FETI<T> &feti, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out)
{
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	math::set(out, T{0});
	for (size_t d = 0; d < L.domain.size(); ++d) {
		for (esint r = 0; r < L.domain[d].B1.nrows; ++r) {
			for (esint c = L.domain[d].B1.rows[r]; c < L.domain[d].B1.rows[r + 1]; ++c) {
				out.vals[L.domain[d].D2C[r]] += L.domain[d].B1.vals[c] * in[d].vals[L.domain[d].B1.cols[c]];
			}
		}
	}
	out.synchronize();
}

template <typename T>
static void insertDomains(FETI<T> &feti, const std::vector<Vector_Dense<T> > &in, Vector_Dual<T> &out)
{
	const typename FETI<T>::EqualityConstraints &L = feti.equalityConstraints;

	math::set(out, T{0});
	for (size_t d = 0; d < L.domain.size(); ++d) {
		for (esint r = 0; r < L.domain[d].B1.nrows; ++r) {
			out.vals[L.domain[d].D2C[r]] += in[d].vals[r];
		}
	}
	out.synchronize();
}

}

#endif /* SRC_FETI_COMMON_APPLYB_H_ */
