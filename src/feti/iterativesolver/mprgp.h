
#ifndef SRC_FETI_ITERATIVESOLVER_MPRGP_H_
#define SRC_FETI_ITERATIVESOLVER_MPRGP_H_

#include "iterativesolver.h"

namespace espreso {

// Modified Proportioning with Reduced Gradient Projection

template <typename T>
class MPRGP: public IterativeSolver<T> {
public:
    MPRGP(FETI<T> &feti);

    void info();
    void solve(const step::Step &step, IterativeSolverInfo &info);

    using IterativeSolver<T>::feti;
    Vector_Dual<T> x, x0, nx, ng, ngg0, xp, g_red, g_free, g_stop;
    Vector_Dual<T> z, p, Fp, g, g0, gg0;
    Vector_Dual<int> free, active;

protected:
    void restrictToFeasibleSet(Vector_Dual<T> &x);
    void updateFreeAndActiveSet(Vector_Dual<int> &free, Vector_Dual<int> &active, const Vector_Dual<T> &x, T epsilon);
    void updateReducedGradient(Vector_Dual<T> &g_red, const Vector_Dual<T> &g, const Vector_Dual<T> &x, double alpha, double epsilon);
    void updateStoppingGradient(Vector_Dual<T> &g_stop, Vector_Dual<T> &g, Vector_Dual<T> &x, double alpha, double epsilon);
    void updateFreeGradient(Vector_Dual<T> &g_free, Vector_Dual<T> &g, Vector_Dual<int> &x, double alpha, double epsilon);
    void multByFree(Vector_Dual<T> &z, Vector_Dual<int> &free);
    double getFeasibleStepLength(Vector_Dual<T> &x, Vector_Dual<T> &p);
    bool stop(const Vector_Dual<T> &x, double epsilon);
};

}


#endif /* SRC_FETI_ITERATIVESOLVER_MPRGP_H_ */
