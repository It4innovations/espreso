/*
 * LinearSolver.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */
//#include <Driver/DissectionSolver.hpp>
#include "LinearSolver.h"
#include "../../basis/utilities/utils.h"



//#include <Eigen/Dense>
//using Eigen::MatrixXd;

using namespace espreso;

LinearSolver::LinearSolver(Instance *instance, const ESPRESOSolver &configuration)
: instance(instance),
  configuration(configuration),
  timeEvalMain("ESPRESO Solver Overal Timing"),
  physics(NULL),
  constraints(NULL),
  cluster(NULL),
  solver(NULL)
{
}


LinearSolver::LinearSolver(const ESPRESOSolver &configuration, OldPhysics &physics, Constraints &constraints)
: instance(NULL),
  configuration(configuration),
  timeEvalMain("ESPRESO Solver Overal Timing"),
  physics(&physics),
  constraints(&constraints)
{
	cluster  = new SuperCluster(configuration, instance);
	solver = new IterSolver(configuration);
}

LinearSolver::~LinearSolver() {

	if (cluster != NULL) {
		delete cluster;
	}

	if (solver != NULL) {
		delete solver;
	}
}


// make full initialization of solver
void LinearSolver::init()
{
	if (cluster != NULL) {
		delete cluster;
	}

	if (solver != NULL) {
		delete solver;
	}

	instance->computeKernels(configuration.regularization, configuration.SC_SIZE);
	if (configuration.method == ESPRESO_METHOD::HYBRID_FETI) {
		instance->assembleB0(configuration.B0_type, instance->N1);
	}
	solver  = new IterSolver	(configuration);
	cluster = new SuperCluster	(configuration, instance);

	init(instance->neighbours);
}

// make partial initialization according to updated matrices
void LinearSolver::update(Matrices matrices)
{
	// TODO update appropriate solver objects and stop steeling matrices! :)

	if (matrices & Matrices::K) {
		// factorization and preconditioners and HFETI preprocessing

		delete cluster;
		delete solver;

		instance->computeKernels(configuration.regularization, configuration.SC_SIZE);
		cluster = new SuperCluster(configuration, instance);
		solver  = new IterSolver(configuration);

		init(instance->neighbours);

	}

	if (matrices & (Matrices::N | Matrices::B1)) { // N is kernel of matrix K
		// updateGGt();
		setup_CreateG_GGt_CompressG();
	}

	if (matrices & (Matrices::B1c | Matrices::f)) {
		// update dual RHS
		//for f:
		// - updated by solve

		//for B1c:
		setup_SetDirichletBoundaryConditions();
	}

	if (matrices & Matrices::B0) {
		// HFETI preprocessing
		if (configuration.method == ESPRESO_METHOD::HYBRID_FETI) {
			instance->assembleB0(configuration.B0_type, instance->N1);
		}
		setup_HTFETI();
	}

	if (matrices & Matrices::B1duplicity) {
		// update duplicity vector
		#pragma omp parallel for
		for (eslocal d = 0; d < cluster->domains.size(); d++) {
			cluster->domains[d]->B1_scale_vec = instance->B1duplicity[d];
		}
	}

	// TODO: remove full re-initialization
	//	delete cluster;
	//	delete solver;
	//	cluster = new Cluster(configuration);
	//	solver = new IterSolver(configuration);
	//	setup();
	//	init(instance->neighbours);
}

// run solver and store primal and dual solution
void LinearSolver::run()
{
	if (
			std::any_of(instance->K.begin(), instance->K.end(), [] (const SparseMatrix &K) { return K.mtype == MatrixType::REAL_UNSYMMETRIC; }) &&
			configuration.solver != ESPRESO_ITERATIVE_SOLVER::GMRES &&
			configuration.solver != ESPRESO_ITERATIVE_SOLVER::BICGSTAB) {

		ESINFO(ERROR) << "Invalid Linear Solver configuration: Only GMRES and BICGSTAB can solve unsymmetric system.";
	}
	Solve(instance->f, instance->primalSolution, instance->dualSolution);
}


void LinearSolver::setup_HTFETI() {
// Setup Hybrid FETI part of the solver

	if (cluster->USE_HFETI == 1) {
		 TimeEvent timeHFETIprec(string("Solver - HFETI preprocessing"));
		 timeHFETIprec.start();

		cluster->SetClusterHFETI();

		 timeHFETIprec.endWithBarrier();
		 timeEvalMain.addEvent(timeHFETIprec);

		ESLOG(MEMORY) << "After HFETI preprocessing process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	}
}

void LinearSolver::setup_LocalSchurComplement() {
// Computation of the Local Schur Complements
//TODO: update for multiple clusters

	if (cluster->USE_KINV == 1) {
		 TimeEvent KSCMem(string("Solver - SC asm. w PARDISO-SC mem [MB]")); KSCMem.startWithoutBarrier(GetProcessMemory_u());
		 TimeEvent timeSolSC2(string("Solver - Schur Complement asm. - using PARDISO-SC"));timeSolSC2.start();
		bool USE_FLOAT = false;
		if (configuration.schur_precision == FLOAT_PRECISION::SINGLE) {
			USE_FLOAT = true;
		}
		cluster->Create_SC_perDomain(USE_FLOAT);
		 timeSolSC2.endWithBarrier(); timeEvalMain.addEvent(timeSolSC2);
		 KSCMem.endWithoutBarrier(GetProcessMemory_u()); //KSCMem.printLastStatMPIPerNode();
		 ESLOG(MEMORY) << "After K inv. process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	} else {
		for (size_t d = 0; d < cluster->domains.size(); d++) {
			cluster->domains[d]->isOnACC = 0;
		}
	}

}

void LinearSolver::createCMat() {
/*

	//int d = 0;


	cluster->SetupPreconditioner();

	SEQ_VECTOR<SparseMatrix> Ml_per_dom (cluster->domains.size());
	SEQ_VECTOR<SparseMatrix> Ml_fin_pd  (cluster->domains.size());

	SEQ_VECTOR<SparseMatrix> Bdec   (cluster->domains.size());
	SEQ_VECTOR<SparseMatrix> Btdec  (cluster->domains.size());



	for (size_t d = 0; d < cluster->domains.size(); d++) {


		cluster->domains[d].B1 = instance->B1[d];

		SparseMatrix B1full = instance->B1[d];
		B1full.ConvertToCSR(0);

		SparseMatrix B1tfull = instance->B1[d];
		B1tfull.ConvertToCSR(0);
		B1tfull.MatTranspose();

		SparseMatrix PrecFull = cluster->domains[d].Prec;

		PrecFull.ConvertDenseToCSR(1);


		SparseMatrix Bt = cluster->domains[d].B1t_DirPr;
		SparseMatrix B;
		Bt.MatTranspose(B);

		Bdec[d] = B;
		Bdec[d].ConvertToCOO(1);

		Btdec[d] = Bt;
		Btdec[d].ConvertToCOO(1);

		SparseMatrix Tmp  = PrecFull;
		Tmp.SetDiagonalOfSymmetricMatrix(0.0);
		Tmp.MatTranspose();
		PrecFull.MatAddInPlace(Tmp,'N',1.0);
		PrecFull.ConvertToCOO(1);
		PrecFull.type='G';

		SEQ_VECTOR <eslocal> prec_map_vec;
		prec_map_vec = cluster->domains[d].B1t_Dir_perm_vec;

		for (size_t i = 0; i < PrecFull.I_row_indices.size(); i++) {
			PrecFull.I_row_indices[i] = prec_map_vec[PrecFull.I_row_indices[i]-1];
			PrecFull.J_col_indices[i] = prec_map_vec[PrecFull.J_col_indices[i]-1];
		}

		PrecFull.rows = cluster->domains[d].B1.cols;
		PrecFull.cols = cluster->domains[d].B1.cols;
		PrecFull.ConvertToCSR(0);

		SparseMatrix TmpBt, Ml;

		TmpBt.MatMat(PrecFull,'N',B1tfull);
		Ml.MatMat(B1full,'N',TmpBt);
		//std::cout << Ml.SpyText();


		//Decompress matrix
		Ml.ConvertToCOO(0);

//		eslocal dl_size = cluster->my_lamdas_indices.size();
//		Ml.rows = dl_size;
//		Ml.cols = dl_size;
//
//
//		for (eslocal i = 0; i < Ml.I_row_indices.size(); i++) {
//			eslocal I 		= Ml.I_row_indices[i] - 1;
//			eslocal J 		= Ml.J_col_indices[i] - 1;
//			eslocal I_dec   = cluster->domains[d].lambda_map_sub_local[I] + 1;
//			eslocal J_dec   = cluster->domains[d].lambda_map_sub_local[J] + 1;
//			Ml.I_row_indices[i] = I_dec;
//			Ml.J_col_indices[i] = J_dec;
//
//		}
//
//		Ml.ConvertToCSR(1);
		Ml_per_dom[d].swap(Ml);
//
//		Bdec[d].rows = dl_size;
//		Btdec[d].cols = dl_size;
//		for (eslocal i = 0; i < Bdec[d].I_row_indices.size(); i++) {
//			eslocal I     = Bdec[d].I_row_indices[i] - 1;
//			eslocal I_dec = cluster->domains[d].lambda_map_sub_local[I] + 1;
//
//			Bdec[d].I_row_indices[i]  = I_dec;
//			Btdec[d].J_col_indices[i] = I_dec;
//		}
//		Bdec[d].ConvertToCSR(1);
//		Btdec[d].ConvertToCSR(1);

	}

	//TODO: can be done as log N
	SparseMatrix Ml_per_cluster;
	Ml_per_cluster = Ml_per_dom[0];
	for (size_t d = 1; d < cluster->domains.size(); d++) {
		Ml_per_cluster.MatAddInPlace(Ml_per_dom[d],'N',1.0);
	}

	SparseMatrix C_per_cluster;

	std::cout << Ml_per_cluster.SpyText();

	for (size_t d = 0; d < cluster->domains.size(); d++) {

		SparseMatrix B1tfull = instance->B1[d];
		SparseMatrix B1full;

		B1tfull.ConvertToCSR(1);
		B1tfull.MatTranspose();
		Esutils::removeDuplicity(B1tfull.CSR_I_row_indices);
		B1tfull.rows = B1tfull.CSR_I_row_indices.size() - 1;

		B1full = B1tfull;
		B1tfull.ConvertToCOO(0);
		B1full.MatTranspose();
		B1full.ConvertToCOO(0);

		SparseMatrix tmp;
		tmp.MatMat(Ml_per_cluster,'N',B1full);
		Ml_fin_pd[d].MatMat(B1tfull,'N',tmp);
		std::cout << Ml_fin_pd[d].SpyText();
		Ml_fin_pd[d].ConvertToCOO(0);

		SparseMatrix Prec = cluster->domains[d].Prec;
		Prec.ConvertDenseToCSR(1);
		SparseMatrix Tmp2 = Prec;
		Tmp2.SetDiagonalOfSymmetricMatrix(0.0);
		Tmp2.MatTranspose();
		Prec.MatAddInPlace(Tmp2,'N',1.0);
		Prec.type='G';

		SparseMatrix A = Prec;
		A.type = 'G';
		SparseMatrix B = Ml_fin_pd[d];
		B.type = 'G';

		A.ConvertCSRToDense(0);
		B.ConvertCSRToDense(0);




//		SparseMatrix Ab = A;
//		SparseMatrix Bb = B;
//
//		int info = 0;
//		info = LAPACKE_dsyev (
//				LAPACK_COL_MAJOR, //int matrix_layout,
//				'N', // char jobz,
//				'U', // char uplo,
//				n, // lapack_int n,
//				&Ab.dense_values[0], //float* a,
//				n, //lapack_int lda,
//				&w[0]); //float* w);
//
//
//		info = 0;
//		info = LAPACKE_dsyev (
//				LAPACK_COL_MAJOR, //int matrix_layout,
//				'N', // char jobz,
//				'U', // char uplo,
//				n, // lapack_int n,
//				&Bb.dense_values[0], //float* a,
//				n, //lapack_int lda,
//				&w[0]); //float* w);
//
//		SparseMatrix Ac = A;
//		SparseMatrix Bc = B;
//
//		info = 0;
//		info = LAPACKE_dsygv (
//				LAPACK_COL_MAJOR, //int matrix_layout,
//				1, //lapack_int itype - itype = 1, the problem type is A*x = lambda*B*x;
//				'V', //char jobz, - If jobz = 'V', then compute eigenvalues and eigenvectors.
//				'U', // char uplo,
//				n, // lapack_int n,
//				&Ac.dense_values[0], // float* a,
//				n, // lapack_int lda,
//				&Bc.dense_values[0], // float* b,
//				n, // lapack_int ldb,
//				&w[0]); //float* w);
//
//		std:fill(w.begin(), w.end(), 0.0);



		int il = 1; // first eigen value to be calculated
		int iv = 4; // last  eigen value to be calculated


		int eig_n = (iv-il+1);
		eslocal n = A.rows;
		SEQ_VECTOR <double> w (eig_n);				// eigen values storage
		SEQ_VECTOR <double> eig_vectors (eig_n*n);	// eigen vectors storage
		SEQ_VECTOR <eslocal> ifail (n);				// dummy
		SEQ_VECTOR <eslocal> m (n);					// dummy
		double tmpd = 0;							// dummy
		double abstol = 0.0;						// dummy



		LAPACKE_dsygvx (
				LAPACK_COL_MAJOR, 	// int matrix_layout,
				1, 					//lapack_int itype,
				'V', 				//char jobz,
				'I', 				//char range,
				'U', 				//char uplo,
				n, 					//lapack_int n,
				&A.dense_values[0], //double* a,
				n, 					//lapack_int lda,
				&B.dense_values[0], //double* b,
				n, 					//lapack_int ldb,
				tmpd, 				//double vl,
				tmpd, 				//double vu,
				il, 				//lapack_int il,
				iv, 				//lapack_int iu,
				abstol, 			//double abstol,
				&m[0],				//lapack_int* m,
				&w[0],				//double* w,
				&eig_vectors[0], 	//double* z,
				n, 					//lapack_int ldz,
				&ifail[0]); 		//lapack_int* ifail);

            int copy_ind_begin = 0;
            for (int i = 0; i < eig_n; i++) {
        		if ( w[i] > 1e-10 ) {
        			copy_ind_begin = n*(i);
        			break;
        		}
        	}

            SparseMatrix Vi;
            Vi.dense_values = std::vector<double> (eig_vectors.begin() + copy_ind_begin, eig_vectors.end() );
            Vi.type = 'G';
            Vi.nnz  = Vi.dense_values.size();
            Vi.rows = n;
			Vi.cols = eig_n - (copy_ind_begin/n);
			Vi.ConvertDenseToCSR(0);

			SparseMatrix BtVi;
			BtVi.MatMat(B1full,'N',Vi);
			BtVi.MatTranspose();

			C_per_cluster.MatAppend(BtVi);

	}

	C_per_cluster.MatTranspose();
	C_per_cluster.ConvertCSRToDense(0);

	SparseMatrix CP_per_cluster;
	SparseMatrix CPt_per_cluster;
	SparseMatrix ACP_per_cluster;

	solver->Projector_l_inv_compG(solver->timeEvalProj, *cluster, C_per_cluster, CP_per_cluster );

	solver->apply_A_l_Mat (solver->timeEvalAppa, *cluster, CP_per_cluster, ACP_per_cluster);

	SparseMatrix Proj_per_cluster;
	SparseMatrix Proj_per_cluster2;

	SparseMatrix Proj_per_cluster_sp;

	CPt_per_cluster = CP_per_cluster;
	CPt_per_cluster.ConvertDenseToCSR(1);
	CPt_per_cluster.MatTranspose();
	CPt_per_cluster.ConvertCSRToDense(1);

	Proj_per_cluster.DenseMatMat(CPt_per_cluster,'N',ACP_per_cluster,'N');

	Proj_per_cluster2.DenseMatMat(CP_per_cluster,'T',ACP_per_cluster,'N'); // TODO : DGEMM nefunguje s transpose turned on - needs to be fixed

	ACP_per_cluster.ConvertDenseToCSR(0);
	CPt_per_cluster.ConvertDenseToCSR(0);

	Proj_per_cluster_sp.MatMat(CPt_per_cluster,'N',ACP_per_cluster);
	Proj_per_cluster_sp.ConvertCSRToDense(0);


	Proj_per_cluster.RemoveLowerDense();
	Proj_per_cluster.mtype = espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

	DenseSolverMKL ProjSolver;

	ProjSolver.ImportMatrix(Proj_per_cluster);

	ProjSolver.Factorization("Geneo projector factorization");




	exit(0);



	SEQ_VECTOR<double> y_out, x_in;

	#pragma omp parallel for
	for (size_t d = 0; d < cluster->domains.size(); d++) {

		SEQ_VECTOR < double > x_in_tmp ( cluster->domains[d].B1_comp_dom.rows, 0.0 );

		for (size_t i = 0; i < cluster->domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_in[ cluster->domains[d].lambda_map_sub_local[i]] * cluster->domains[d].B1_scale_vec[i]; // includes B1 scaling

		cluster->domains[d].B1t_DirPr.MatVec (x_in_tmp, cluster->x_prim_cluster1[d], 'N');
		//cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
		cluster->domains[d].Prec.DenseMatVec(cluster->x_prim_cluster1[d], cluster->x_prim_cluster2[d],'N');



	}



	std::fill( cluster->compressed_tmp.begin(), cluster->compressed_tmp.end(), 0.0);
	SEQ_VECTOR < double > y_out_tmp;
	for (size_t d = 0; d < cluster->domains.size(); d++) {
		y_out_tmp.resize( cluster->domains[d].B1_comp_dom.rows );


		cluster->domains[d].B1t_DirPr.MatVec (cluster->x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

	}

	//solver->All_Reduce_lambdas_compB(cluster, cluster->compressed_tmp, y_out);

*/
}

void LinearSolver::setup_Preconditioner() {
// Load Matrix K, Regularization
		 TimeEvent timeRegKproc(string("Solver - Setup preconditioners")); timeRegKproc.start();
		 ESLOG(MEMORY) << "Before - Setup preconditioners - process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 TimeEvent KregMem(string("Solver - Setup preconditioners mem. [MB]")); KregMem.startWithoutBarrier(GetProcessMemory_u());
		 ESINFO(PROGRESS3) << "Setup preconditioners";

		cluster->SetupPreconditioner();

		 KregMem.endWithoutBarrier(GetProcessMemory_u()); //KregMem.printLastStatMPIPerNode();

		 ESLOG(MEMORY) << "After - Setup preconditioners " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
		 timeRegKproc.endWithBarrier();
		 timeEvalMain.addEvent(timeRegKproc);
}

void LinearSolver::setup_FactorizationOfStiffnessMatrices() {
// K Factorization
		 TimeEvent timeSolKproc(string("Solver - K factorization")); timeSolKproc.start();
		 TimeEvent KFactMem(string("Solver - K factorization mem. [MB]")); KFactMem.startWithoutBarrier(GetProcessMemory_u());
		 ESINFO(PROGRESS3) << "Factorize K";

		cluster->SetupKsolvers();

		 KFactMem.endWithoutBarrier(GetProcessMemory_u()); //KFactMem.printLastStatMPIPerNode();
		 ESLOG(MEMORY) << "After K solver setup process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
		 timeSolKproc.endWithBarrier();
		 timeEvalMain.addEvent(timeSolKproc);
}


void LinearSolver::setup_SetDirichletBoundaryConditions() {
// Set Dirichlet Boundary Condition

		 TimeEvent timeSetInitialCondition(string("Solver - Set Dirichlet Boundary Condition"));
		 timeSetInitialCondition.start();

		#pragma omp parallel for
		for (int d = 0; d < cluster->domains.size(); d++) {
			cluster->domains[d]->vec_c  = instance->B1c[d];
			cluster->domains[d]->vec_lb = instance->LB[d];
		}

		 timeSetInitialCondition.endWithBarrier();
		 timeEvalMain.addEvent(timeSetInitialCondition);
}

void LinearSolver::setup_CreateG_GGt_CompressG() {

		 TimeEvent timeSolPrec(string("Solver - FETI Preprocessing")); timeSolPrec.start();

		 ESLOG(MEMORY) << "Solver Preprocessing - HFETI with regularization from K matrix";
		 ESLOG(MEMORY) << "process " << environment->MPIrank << " uses "<< Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 TimeEvent G1_perCluster_time("Setup G1 per Cluster time   - preprocessing"); G1_perCluster_time.start();

		 TimeEvent G1_perCluster_mem("Setup G1 per Cluster memory - preprocessing"); G1_perCluster_mem.startWithoutBarrier(GetProcessMemory_u());
		cluster->Create_G_perCluster();

		 G1_perCluster_time.end(); G1_perCluster_time.printStatMPI();
		 G1_perCluster_mem.endWithoutBarrier(GetProcessMemory_u()); G1_perCluster_mem.printStatMPI();

		 ESLOG(MEMORY) << "Created G1 per cluster";
		 ESLOG(MEMORY) << "Before HFETI create GGt process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 TimeEvent solver_Preprocessing_time("Setup GGt time   - preprocessing"); solver_Preprocessing_time.start();
		 TimeEvent solver_Preprocessing_mem("Setup GGt memory - preprocessing"); solver_Preprocessing_mem.start();

		solver->Preprocessing(*cluster);

		 solver_Preprocessing_time.end(); solver_Preprocessing_time.printStatMPI();
		 solver_Preprocessing_mem.end();  solver_Preprocessing_mem.printStatMPI();

		 ESLOG(MEMORY) << "Create GGtInv";
		 ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 TimeEvent solver_G1comp_time("Setup G1 compression time   - preprocessing"); solver_G1comp_time.start();
		 TimeEvent solver_G1comp_mem("Setup G1 compression memory - preprocessing");  solver_G1comp_mem.start();

		cluster->Compress_G1(); // Compression of Matrix G1 to work with compressed lambda vectors

		 solver_G1comp_time.end(); solver_G1comp_time.printStatMPI();
		 solver_G1comp_mem.end();  solver_G1comp_mem.printStatMPI();
		 ESLOG(MEMORY) << "G1 compression";
		 ESLOG(MEMORY) << "process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
		 ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 timeSolPrec.endWithBarrier(); timeEvalMain.addEvent(timeSolPrec);


}



void LinearSolver::setup_InitClusterAndSolver( )
{

	cluster->SetupCommunicationLayer();


	// *** Iter Solver Set-up ***************************************************************************************************************************
	solver->CG_max_iter 		= configuration.iterations;
	solver->USE_GGtINV  		= 1;
	solver->epsilon 			= configuration.epsilon;
	solver->USE_PREC 			= configuration.preconditioner;

	solver->USE_KINV 			= configuration.use_schur_complement ? 1 : 0;
	solver->PAR_NUM_THREADS 	= environment->PAR_NUM_THREADS;
	solver->SOLVER_NUM_THREADS  = environment->SOLVER_NUM_THREADS;

	switch (configuration.method) {
	case ESPRESO_METHOD::TOTAL_FETI:
		solver->USE_HFETI = false;
		break;
	case ESPRESO_METHOD::HYBRID_FETI:
		solver->USE_HFETI = true;
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unsupported FETI METHOD";
	}


}

// TODO: const parameters
void LinearSolver::init(const std::vector<int> &neighbours)
{
	if (physics != NULL) {
		instance = new Instance(physics->K.size(), neighbours);
		instance->computeKernelsCallback = physics->computeKernelsCallback;
		instance->assembleB0Callback = physics->assembleB0Callback;
		instance->clustersMap = constraints->continuityMap;

		physics->K.swap(instance->K);
		for (size_t i = 0; i < physics->K.size(); i++) {
			instance->K[i].mtype = physics->mtype;
		}

		instance->computeKernelsCallback = [&] (REGULARIZATION regularization, size_t scSize) {
			physics->K.swap(instance->K);
			physics->computeKernelsCallback(regularization, scSize);
			physics->K.swap(instance->K);
			physics->R1.swap(instance->N1);
			physics->R2.swap(instance->N2);
			physics->RegMat.swap(instance->RegMat);
		};

		instance->computeKernels(configuration.regularization, configuration.SC_SIZE);

		constraints->B1.swap(instance->B1);
		constraints->B1c.swap(instance->B1c);
		constraints->B1subdomainsMap.swap(instance->B1subdomainsMap);
		constraints->B1clustersMap.swap(instance->B1clustersMap);
		constraints->B1duplicity.swap(instance->B1duplicity);

		if (configuration.method == ESPRESO_METHOD::HYBRID_FETI) {
			instance->assembleB0Callback = [&] (B0_TYPE type, const std::vector<SparseMatrix> &kernels) {
				physics->assembleB0Callback(type, kernels);
				constraints->B0.swap(instance->B0);
				constraints->B0subdomainsMap.swap(instance->B0subdomainsMap);
			};

			instance->assembleB0(configuration.B0_type, instance->N1);
		}

		constraints->LB.swap(instance->LB);
		constraints->inequality.swap(instance->inequality);
		constraints->inequalityC.swap(instance->inequalityC);

		cluster->instance  = instance;
		cluster->init();
	}


	//mkl_cbwr_set(MKL_CBWR_COMPATIBLE);

	// Overall Linear Solver Time measurement structure
	timeEvalMain.totalTime.startWithBarrier();



	// *** Initialize the Cluster and Solver structures
	setup_InitClusterAndSolver();




	// *** Setup B0 matrix
	// setup_B0Matrices();

	// *** Setup B1 matrix
	// setup_B1Matrices();

	// *** Setup R matrix
	// setup_KernelMatrices();

	// *** Set Dirichlet Boundary Condition
	// setup_SetDirichletBoundaryConditions();

	// *** if TFETI is used or if HTFETI and analytical kernel are used we can compute GGt here - between solution in terms of peak memory
	if (  !(cluster->USE_HFETI == 1 && configuration.regularization == REGULARIZATION::NULL_PIVOTS)  ) {
		setup_CreateG_GGt_CompressG();
	}

	// *** setup all preconditioners
	setup_Preconditioner();

	// *** Computation of the Schur Complement
	setup_LocalSchurComplement();

	// *** K Factorization
	setup_FactorizationOfStiffnessMatrices();

	// *** Setup Hybrid FETI part of the solver
	setup_HTFETI();

	// *** For algebraic kernels, GGt needs to be computed after HTFETI preprocessing
	if (cluster->USE_HFETI == 1 && configuration.regularization == REGULARIZATION::NULL_PIVOTS) {
		setup_CreateG_GGt_CompressG();
	}

	// Setup Conj Projector

	if ( configuration.conj_projector == CONJ_PROJECTOR::GENEO ) {
		// solver->CreateConjProjector(*cluster);
	}

	// Cleanup of unnecessary objects
	// TODO: This can be a problem in some cases - need to be verified
	cluster->my_lamdas_map_indices.clear();

	#pragma omp parallel for
	for (size_t d = 0; d < cluster->domains.size(); d++) {
		cluster->domains[d]->B1.Clear();
	}

	ESLOG(MEMORY) << "End of preprocessing - process " << environment->MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

}

void LinearSolver::Solve( std::vector < std::vector < double > >  & f_vec,
		                  std::vector < std::vector < double > >  & prim_solution)
{

	std::vector < std::vector < double > > dual_solution;
	Solve(f_vec, prim_solution, dual_solution);
}

void LinearSolver::Solve( std::vector < std::vector < double > >  & f_vec,
		                  std::vector < std::vector < double > >  & prim_solution,
		                  std::vector < std::vector < double > > & dual_solution) {

	 	 TimeEvent timeSolCG(string("Solver - CG Solver runtime"));
	 	 timeSolCG.start();

	 solver->Solve    ( *cluster, f_vec, prim_solution, dual_solution );

	 	 timeSolCG.endWithBarrier();
	 	 timeEvalMain.addEvent(timeSolCG);

}

void LinearSolver::Postprocessing( ) {

}

void LinearSolver::finilize() {

	// Show Linear Solver Runtime Evaluation
	solver->preproc_timing.printStatsMPI();
	solver->timing.printStatsMPI();
	solver->postproc_timing.printStatsMPI();
	solver->timeEvalAppa.printStatsMPI();
	solver->timeEvalProj.printStatsMPI();

	if ( solver->USE_PREC != ESPRESO_PRECONDITIONER::NONE ) solver->timeEvalPrec.printStatsMPI();

	//TODO: Fix timing:  if ( cluster->USE_HFETI == 1 ) cluster->ShowTiming();

	 timeEvalMain.totalTime.endWithBarrier();
	 timeEvalMain.printStatsMPI();
}

void LinearSolver::CheckSolution( vector < vector < double > > & prim_solution ) {
    // *** Solutin correctnes test **********************************************************************************************

	//	double max_v = 0.0;
//		for (eslocal i = 0; i < number_of_subdomains_per_cluster; i++)
//			for (size_t j = 0; j < prim_solution[i].size(); j++)
//				if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );
//
//	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.startWithoutBarrier(0.0); max_sol_ev.endWithoutBarrier(max_v);
//
//	double max_vg;
//	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, environment->MPICommunicator );
//	ESINFO(DETAILS) << "Maxvalue in solution = " << std::setprecision(12) << max_vg;

	//max_sol_ev.printLastStatMPIPerNode(max_vg);
	// *** END - Solutin correctnes test ******************************************************************************************

}

void LinearSolver::Preprocessing( std::vector < std::vector < eslocal > > & lambda_map_sub) {

}
