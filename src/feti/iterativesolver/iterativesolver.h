
#ifndef SRC_FETI_ITERATIVESOLVER_ITERATIVESOLVER_H_
#define SRC_FETI_ITERATIVESOLVER_ITERATIVESOLVER_H_

#include "feti/feti.h"
#include "feti/common/vector_dual.h"
#include "feti/projector/vector_kernel.h"

namespace espreso {

struct IterativeSolverInfo {
    enum class ERROR: int {
        OK = 0,
        STAGNATION,
        MAX_ITERATIONS_REACHED,
        INVALID_DATA,
        CONVERGENCE_ERROR,
    };

    ERROR error = ERROR::OK;
    size_t iterations = 0;
    bool converged = false;

    struct Norm {
        struct Dual {
            double absolute, relative, arioli, initial, ksi, criteria;
        } dual;
        double primal;
    } norm;

    struct Time {
        double current, total;
    } time;

    struct Stagnation {
        std::vector<double> buffer;
        int p = 0;
    } stagnation;
};

template <typename T>
class IterativeSolver {
public:
    static IterativeSolver<T>* create(FETI<T> &feti, const step::Step &step);

    template <typename Vector>
    void resize(Vector &v)
    {
        v.resize();
    }

    template <typename Vector, typename... Other>
    void resize(Vector &v, Other&... other)
    {
        v.resize();
        resize(other...);
    }

    IterativeSolver(FETI<T> &feti): feti(feti)
    {
        iKfBtL.resize(feti.K.size());
        Ra.resize(feti.K.size());

        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            iKfBtL[d].resize(feti.K[d].nrows);
            Ra[d].resize(feti.K[d].nrows);
        }
    }

    virtual ~IterativeSolver() {}

    virtual void info() =0;
    virtual void set(const step::Step &step) { }
    virtual void update(const step::Step &step) =0;
    virtual void solve(const step::Step &step, IterativeSolverInfo &info) =0;

    void setInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const T &ww);
    void updateInfo(IterativeSolverInfo &info, const FETIConfiguration &configuration, const T &ww, const T &psi, const T &ry);
    void reconstructSolution(const Vector_Dual<T> &l, const Vector_Dual<T> &r, const step::Step &step);
    void reconstructSolution(const Vector_Dual<T> &l, const Vector_Kernel<T> &rbm, const step::Step &step);

    FETI<T> &feti;

    std::vector<Vector_Dense<T> > iKfBtL, Ra;

protected:
    void print(const step::Step &step);
};

}

#endif /* SRC_FETI_ITERATIVESOLVER_ITERATIVESOLVER_H_ */

