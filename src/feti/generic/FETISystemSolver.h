/*
 * LinearSolver.h
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#ifndef SOLVER_GENERIC_LINEARSOLVER_H_
#define SOLVER_GENERIC_LINEARSOLVER_H_

#include <vector>

#include "physics/system/linearsystem.h"

namespace espreso {

class FETISolverData;
struct FETIDataHolder;
struct FETIConfiguration;
class MatrixCSRFETI;
class MatrixIJVFETI;
class MatrixDenseFETI;
class VectorDenseFETI;
class VectorsDenseFETI;
class TimeEval;

class FETISystemSolver: public SystemSolver {
public:
	FETISystemSolver(FETIConfiguration &configuration, FETISolverData &data);

	void init();
	void update();
	void solve();

	double& precision();

	void insertK(FETIConfiguration &configuration, const MatrixCSRFETI &K, const MatrixCSRFETI &origK, const MatrixDenseFETI &N1, const MatrixDenseFETI &N2, const MatrixCSRFETI &RegMat);
	void insertB1(const MatrixIJVFETI &B1Dirichlet, const VectorDenseFETI &c, const MatrixIJVFETI &B1Gluing, const VectorDenseFETI &duplication, const MatrixIJVFETI &B1Inequality, const VectorDenseFETI &gap, const std::vector<esint> &B1Map);
	void insertB0(const MatrixIJVFETI &B0);
	void insertRHS(const VectorsDenseFETI &f);

	void update(FETIConfiguration &configuration);
	void solve(FETIConfiguration &configuration, VectorsDenseFETI &x, VectorsDenseFETI &y);

	virtual ~FETISystemSolver();

//	void setup();

	void init(const std::vector<int> &neighbors, FETIConfiguration &configuration);

	void Preprocessing( std::vector < std::vector < esint > > & lambda_map_sub );

	int Solve( std::vector < std::vector <double > >	& f_vec, std::vector < std::vector < double > > & prim_solution );
	int Solve( std::vector < std::vector <double > >	& f_vec, std::vector < std::vector < double > > & prim_solution, std::vector < std::vector < double > > & dual_solution );

	void Postprocessing ();

	void finalize();

	void CheckSolution( std::vector < std::vector < double > > & prim_solution );

	FETIConfiguration &configuration;
	TimeEval *timeEvalMain;

private:

	FETISolverData &_data;

	FETIDataHolder *_inner;
	void setup_HTFETI();

	void setup_LocalSchurComplement(FETIConfiguration &configuration);
	void setup_Preconditioner();
	void setup_FactorizationOfStiffnessMatrices();
	void setup_SetDirichletBoundaryConditions();

	void setup_CreateG_GGt_CompressG();
	void setup_InitClusterAndSolver(FETIConfiguration &configuration);
};

}

#endif /* SOLVER_GENERIC_LINEARSOLVER_H_ */
