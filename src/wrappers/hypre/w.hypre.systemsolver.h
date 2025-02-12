
#ifndef SRC_WRAPPERS_HYPRE_W_HYPRE_SYSTEMSOLVER_H_
#define SRC_WRAPPERS_HYPRE_W_HYPRE_SYSTEMSOLVER_H_

#include "physics/system/linearsystem.h"

namespace espreso {

struct HYPREConfiguration;
struct HYPRESolverData;
struct HYPREDataHolder;

class HYPRESystemSolver: public SystemSolver {
public:
    HYPRESystemSolver(HYPREConfiguration &configuration, HYPRESolverData &data);

    void init();
    void update();
    void solve();

    double& precision();

    ~HYPRESystemSolver();

    HYPREConfiguration &configuration;

protected:
    esint _roffset, _nrows;

    HYPRESolverData &_data;

    HYPREDataHolder *_inner;
};

}

#endif /* SRC_WRAPPERS_HYPRE_W_HYPRE_SYSTEMSOLVER_H_ */
