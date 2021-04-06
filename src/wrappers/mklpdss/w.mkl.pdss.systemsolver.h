
#ifndef SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_SYSTEMSOLVER_H_
#define SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_SYSTEMSOLVER_H_

#include "physics/system/linearsystem.h"

namespace espreso {

struct MKLPDSSSolverData;
struct MKLPDSSConfiguration;
struct MKLPDSSDataHolder;

class MKLPDSSSystemSolver: public SystemSolver {
public:
	MKLPDSSSystemSolver(MKLPDSSConfiguration &configuration, MKLPDSSSolverData &data);

	void init();
	void update();
	void solve();

	double& precision() { return _precision; }

	~MKLPDSSSystemSolver();

	MKLPDSSConfiguration &configuration;

protected:
	void call(esint phase);

	esint _roffset, _nrows;
	double _precision; // dummy

	MKLPDSSSolverData &_data;

	MKLPDSSDataHolder *_inner;
};

}


#endif /* SRC_WRAPPERS_MKLPDSS_W_MKL_PDSS_SYSTEMSOLVER_H_ */
