
#ifndef SRC_FETI_PROJECTOR_PROJECTOR_H_
#define SRC_FETI_PROJECTOR_PROJECTOR_H_

#include "feti/feti.h"
#include "feti/common/vector_dual.h"
#include "feti/projector/vector_kernel.h"

namespace espreso {

template <typename T>
struct Projector {
    static Projector<T>* create(FETI<T> &feti, const step::Step &step);

    Projector(FETI<T> &feti);
    virtual ~Projector() {}

    virtual void set(const step::Step &step) =0;
    virtual void update(const step::Step &step) =0;

    void info();

    // I - Q
    void apply(const Vector_Dual<T> &x, Vector_Dual<T> &y);
    void apply_GtintGGt(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
    void apply_R(const Vector_Kernel<T> &x, std::vector<Vector_Dense<T> > &y);
    void apply_RinvGGtG(const Vector_Dual<T> &x, std::vector<Vector_Dense<T> > &y);
    void apply_invU(const Vector_Kernel<T> &x, Vector_Kernel<T> &y);
    void apply_invL(const Vector_Kernel<T> &x, Vector_Kernel<T> &y);
    void apply_GtinvU(const Vector_Kernel<T> &x, Vector_Dual<T> &y);
    void apply_invLG(const Vector_Dual<T> &x, Vector_Kernel<T> &y);

    Vector_Kernel<T> e;

protected:
    enum class GGT_TYPE { GGT, LU, NONE };

    struct Kernel {
        static int roffset, rsize, total;
        int offset, size;
    };

    void _apply_G(const Vector_Dual<T> &in, Vector_Kernel<T> &out);
    void _apply_invGGt(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
    void _apply_invL(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
    void _apply_invU(const Vector_Kernel<T> &in, Vector_Dense<T> &out);
    void _apply_Gt(const Vector_Dense<T> &in, const T &alpha, Vector_Dual<T> &out);
    void _apply_R(const Vector_Dense<T> &in, std::vector<Vector_Dense<T> > &out);

    void _print(const step::Step &step);

    FETI<T> &feti;

    GGT_TYPE GGTtype;
    std::vector<Kernel> kernel;

    Matrix_CSR<T> G, Gt; // G has just pointers to Gt
    Matrix_CSR<T> GGt;   // global object on all CPUs
    Matrix_Dense<T> invGGt, invL, invU; // GGt or U, L

    Vector_Kernel<T> Gx; // we need whole vector
    Vector_Dense<T> iGGtGx; // only local part is sufficient
};

}

#endif /* SRC_FETI_PROJECTOR_PROJECTOR_H_ */
