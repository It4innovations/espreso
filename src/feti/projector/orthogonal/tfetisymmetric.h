
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
	OrthogonalTFETISymmetric(FETI<T> *feti);
	~OrthogonalTFETISymmetric();

	void info();
	void update();

	void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
	void applyGtInvGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
	void applyRInvGGtG(const Vector_Dual<T> &x, Vector_FETI<Vector_Dense, T> &y);

protected:
	void _setG();
	void _setSparseGGt();
	void _setDenseGGt();
	void _updateG();
	void _updateSparseGGt();
	void _updateDenseGGt();

	void _applyG(const Vector_Dual<T> &in, Vector_Kernel<T> &out);
	void _applyInvGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
	void _applyGt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out);
	void _applyR(const Vector_Dense<T> &in, Vector_FETI<Vector_Dense, T> &out);

	void _print();

	Matrix_CSR<T> G, GGt;
	Matrix_Dense<T> invGGt;

	Vector_Kernel<T> Gx; // we need whole vector
	Vector_Dense<T> iGGtGx; // only local part is sufficient

	size_t GGtOffset, GGtSize;
	std::vector<esint> Roffset;
	std::vector<std::pair<esint, esint> > Goffset; // offset to G for each LMAL
	std::vector<std::vector<std::pair<esint, esint> > > nKernels; // n, offset
	std::vector<std::vector<T> > sBuffer, rBuffer;
	DirectSolver<T, Matrix_CSR> GGtSolver;
};

}




#endif /* SRC_FETI_PROJECTOR_ORTHOGONAL_TFETISYMMETRIC_H_ */
