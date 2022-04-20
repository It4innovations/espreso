
#ifndef SRC_AXFETI_PROJECTOR_ORTHOGONAL_TFETISYMMETRIC_H_
#define SRC_AXFETI_PROJECTOR_ORTHOGONAL_TFETISYMMETRIC_H_

#include "axfeti/projector/projector.h"

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
class OrthogonalTFETISymmetric: public Projector<T> {
public:
	OrthogonalTFETISymmetric(AX_FETI<T> *feti);
	~OrthogonalTFETISymmetric();

	void info();
	void update();

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
	void applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
	void applyInvGGtG(const Vector_Dual<T> &x, Vector_Kernel<T> &y);

	Matrix_CSR<T> G, GGt;
	Matrix_Dense<T> invGGt;

	Vector_Kernel<T> Gx; // we need whole vector
	Vector_Dense<T> iGGtGx; // only local part is sufficient

	size_t GGtOffset, GGtSize;
	std::vector<esint> Roffset;
	std::vector<std::pair<esint, esint> > Goffset; // offset to G for each LMAL
	std::vector<std::vector<std::pair<esint, esint> > > nKernels; // n, offset
	std::vector<std::vector<T> > sBuffer, rBuffer;
};

}




#endif /* SRC_AXFETI_PROJECTOR_ORTHOGONAL_TFETISYMMETRIC_H_ */
