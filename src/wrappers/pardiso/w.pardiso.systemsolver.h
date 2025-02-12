
#ifndef SRC_WRAPPERS_PARDISO_W_PARDISO_SYSTEMSOLVER_H_
#define SRC_WRAPPERS_PARDISO_W_PARDISO_SYSTEMSOLVER_H_

#include "physics/system/linearsystem.h"

namespace espreso {

struct PARDISOSolverData;
struct PARDISOConfiguration;
class VectorDense;
class MatrixCSR;


class PARDISOSystemSolver: public SystemSolver {
public:
    PARDISOSystemSolver(PARDISOConfiguration &configuration, PARDISOSolverData &data);

    void init();
    void update();
    void solve();

    double& precision() { return _precision; }
    ~PARDISOSystemSolver();

    PARDISOConfiguration &configuration;

protected:
    void call(esint phase);

    double _precision; // dummy

    PARDISOSolverData &_data;

    MatrixCSR *_K;
    VectorDense *_f, *_x;
};

}



#endif /* SRC_WRAPPERS_PARDISO_W_PARDISO_SYSTEMSOLVER_H_ */
