
#ifndef SRC_WRAPPERS_WSMP_W_WSMP_SYSTEMSOLVER_H_
#define SRC_WRAPPERS_WSMP_W_WSMP_SYSTEMSOLVER_H_

#include "physics/system/linearsystem.h"

namespace espreso {

struct WSMPSolverData;
struct WSMPConfiguration;
struct WSMPDataHolder;

class WSMPSystemSolver: public SystemSolver {
public:
    WSMPSystemSolver(WSMPConfiguration &configuration, WSMPSolverData &data);

    void init();
    void update();
    void solve();

    double& precision() { return _precision; }

    ~WSMPSystemSolver();

    WSMPConfiguration &configuration;

protected:
    esint _roffset, _nrows;
    double _precision; // dummy

    WSMPSolverData &_data;

    WSMPDataHolder *_inner;
};

}

#endif /* SRC_WRAPPERS_WSMP_W_WSMP_SYSTEMSOLVER_H_ */
