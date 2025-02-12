#ifndef SRC_WRAPPERS_SUPERLU_W_SUPERLU_SYSTEMSOLVER_H_
#define SRC_WRAPPERS_SUPERLU_W_SUPERLU_SYSTEMSOLVER_H_

#include "physics/system/linearsystem.h"

namespace espreso {

struct SuperLUSolverData;
struct SuperLUConfiguration;
struct SuperLUDataHolder;

class SuperLUSystemSolver: public SystemSolver {
public:
    SuperLUSystemSolver(SuperLUConfiguration &configuration, SuperLUSolverData &data);

    void init();
    void update();
    void solve();

    double& precision() { return _precision; }

    ~SuperLUSystemSolver();

    SuperLUConfiguration &configuration;

protected:
    esint _roffset, _nrows;
    double _precision; // dummy

    SuperLUSolverData &_data;

    SuperLUDataHolder *_inner;
};

}

#endif /* SRC_WRAPPERS_SUPERLU_W_SUPERLU_SYSTEMSOLVER_H_ */
