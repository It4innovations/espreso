
#ifndef SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_
#define SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_

#include "linearsystem.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "feti/feti.h"
#include "math/physics/matrix_feti.h"
#include "math/physics/matrix_distributed.h"
#include "math/physics/vector_feti.h"

namespace espreso {

template <typename T, class Physics>
struct FETILinearSystemSolver: LinearSystemSolver<T> {

	FETILinearSystemSolver(FETIConfiguration &configuration)
	: feti(configuration)
	{
		LinearSystemSolver<T>::A = &feti.K;
		LinearSystemSolver<T>::x = &feti.x;
		LinearSystemSolver<T>::b = &feti.f;
		LinearSystemSolver<T>::dirichlet = &dirichlet;
	}

	~FETILinearSystemSolver()
	{

	}

	void set(step::Step &step)
	{
		eslog::startln("FETI: SETTING LINEAR SYSTEM", "FETI[SET]");
		setEqualityConstraints(step);
		eslog::checkpointln("FETI: SET B1");
		setKernel(step);
		eslog::checkpointln("FETI: SET KERNELS");
		eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
		feti.set(step);
		eslog::endln("FETI: LINEAR SYSTEM SET");
	}

	void update(step::Step &step)
	{
		eslog::startln("FETI: UPDATING LINEAR SYSTEM", "FETI[UPDATE]");
		updateEqualityConstraints(step);
		eslog::checkpointln("FETI: UPDATE B1");
		updateKernel(step);
		eslog::checkpointln("FETI: UPDATE KERNELS");
		feti.update(step);
		if (info::ecf->output.print_matrices) {
			eslog::storedata(" STORE: system/{K, f}\n");
			math::store(feti.K, utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
			math::store(feti.f, utils::filename(utils::debugDirectory(step) + "/system", "f").c_str());
		}
		eslog::endln("FETI: LINEAR SYSTEM UPDATED");
	}

	bool solve(step::Step &step)
	{
		eslog::startln("FETI: RUN LINEAR SYSTEM", "FETI[SOLVE]");
		if (feti.solve(step)) {
			if (info::ecf->output.print_matrices) {
				eslog::storedata(" STORE: system/{x}\n");
				math::store(feti.x, utils::filename(utils::debugDirectory(step) + "/system", "x").c_str());
			}
			eslog::endln("FETI: LINEAR SYSTEM SOLVED");
			return true;
		}
		eslog::endln("FETI: LINEAR SYSTEM SOLVED");
		return false;
	}

private:
	void setEqualityConstraints(step::Step &step);
	void updateEqualityConstraints(step::Step &step);

	void setKernel(step::Step &step);
	void updateKernel(step::Step &step);

	Vector_Distributed<Vector_Sparse, T> dirichlet;

	FETI<T> feti;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_FETISOLVER_H_ */
