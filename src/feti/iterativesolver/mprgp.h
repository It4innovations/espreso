
#ifndef SRC_FETI_ITERATIVESOLVER_MPRGP_H_
#define SRC_FETI_ITERATIVESOLVER_MPRGP_H_

#include "iterativesolver.h"

namespace espreso {

// Modified Proportioning with Reduced Gradient Projection

struct MPRGPSolverInfo: IterativeSolverInfo {
    bool print;
    int n_cg, n_mixed, n_gproj, n_hess;
};

template <typename T>
class MPRGP: public IterativeSolver<T> {
public:
    MPRGP(FETI<T> &feti);

    void info();

    void update(const step::Step &step)
    {
        IterativeSolver<T>::resize(b, x, x0, nx, ng, ngg0, xp, g_red, g_free, g_stop);
        IterativeSolver<T>::resize(z, p, Fp, g, g0, gg0);
        IterativeSolver<T>::resize(free, active);
        IterativeSolver<T>::resize(Gx);
    }

    void solve(const step::Step &step, IterativeSolverInfo &info);

    void run(const step::Step &step, MPRGPSolverInfo &info, double alpha, std::function<void(Vector_Dual<T> &in, Vector_Dual<T> &out)> H, std::function<bool(const Vector_Dual<T> &x, const Vector_Dual<T> &g_stop)> stop);

    using IterativeSolver<T>::feti;
    Vector_Dual<T> b, x, x0, nx, ng, ngg0, xp, g_red, g_free, g_stop;
    Vector_Dual<T> z, p, Fp, g, g0, gg0;
    Vector_Dual<int> free, active;
    Vector_Kernel<T> Gx;

    void restrictToFeasibleSet(Vector_Dual<T> &x);
    void updateFreeAndActiveSet(Vector_Dual<int> &free, Vector_Dual<int> &active, const Vector_Dual<T> &x);
    void updateReducedGradient(Vector_Dual<T> &g_red, const Vector_Dual<T> &g, const Vector_Dual<T> &x, double alpha);
    void updateStoppingGradient(Vector_Dual<T> &g_stop, Vector_Dual<T> &g, Vector_Dual<T> &x, double alpha);
    void updateFreeGradient(Vector_Dual<T> &g_free, Vector_Dual<T> &g, Vector_Dual<int> &x);
    void multByFree(Vector_Dual<T> &z, Vector_Dual<int> &free);
    double getFeasibleStepLength(Vector_Dual<T> &x, Vector_Dual<T> &p);
};

}


#endif /* SRC_FETI_ITERATIVESOLVER_MPRGP_H_ */
