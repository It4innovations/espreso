
#ifndef SRC_FETI_PROJECTOR_ORTHOGONAL_TFETISYMMETRIC_H_
#define SRC_FETI_PROJECTOR_ORTHOGONAL_TFETISYMMETRIC_H_

#include "feti/projector/projector.h"

#include <unordered_map>

namespace espreso {

/*
 * R: KxR : block diagonal
 * B: LxK : from primal to dual
 *
 * e = Rt * f : R
 * G = Rt * Bt: RxL: from dual to kernels
 *
 * y = Q * x = Gt * inv(GGt) * G * x
 *
 *     Gx = G * x         :: x      -> Gx     :      L -> R
 * iGGtGx = inv(GGt) * Gx :: Gx     -> iGGtGx : totalR -> R
 *      y = Gt * iGGtGx   :: iGGtGx -> y      :      R -> L
 */

template <typename T>
struct OrthogonalTFETISymmetric: public Projector<T> {
	OrthogonalTFETISymmetric(FETI<T> &feti);
	~OrthogonalTFETISymmetric();

	void info();
	void update(const step::Step &step);

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
	void applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
	void applyRInvGGtG(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y);

protected:
	void _computeDualGraph();
	void _setG();
	void _setGGt();
	void _updateG();
	void _updateGGt();

	void _applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out);
	void _applyInvGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
	void _applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out);
	void _applyR(const Vector_Dense<T> &in, Vector_FETI<Vector_Dense, T> &out);

	void _print(const step::Step &step);

	using Projector<T>::feti;
	using Projector<T>::e;

	Matrix_CSR<T> G, Gt, GGt;
	Matrix_Dense<T> invGGt;

	Vector_Kernel<T> Gx; // we need whole vector
	Vector_Dense<T> iGGtGx; // only local part is sufficient

	size_t GGtOffset, GGtSize, GGtDataOffset, GGtNnz;

	std::vector<std::vector<esint> > dualGraph;

//	struct DomainInfo { esint domain, localOffset, globalOffset, kernels, size; };
//	std::unordered_map<esint, DomainInfo> domainInfo;
};

}




#endif /* SRC_FETI_PROJECTOR_ORTHOGONAL_TFETISYMMETRIC_H_ */
