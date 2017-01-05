/*
 * LinearSolver.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#include "LinearSolver.h"

//#include <Eigen/Dense>
//using Eigen::MatrixXd;

using namespace espreso;

LinearSolver::~LinearSolver() {
	// TODO Auto-generated destructor stub
}

void LinearSolver::setup() {

	SINGULAR 	= physics.singular();

	if (!config::solver::KEEP_FACTORS)
		KEEP_FACTORS = false; // only suported by MKL Pardiso so far
	else
		KEEP_FACTORS = true;

    // ***************************************************************************************************************************
	// Cluster structure  setup
	if ( SINGULAR )
		cluster.USE_DYNAMIC		= 0;
	else
		cluster.USE_DYNAMIC		= 1;

	switch (config::solver::FETI_METHOD) {
	case config::solver::FETI_METHODalternative::TOTAL_FETI:
		cluster.USE_HFETI = false;
		break;
	case config::solver::FETI_METHODalternative::HYBRID_FETI:
		cluster.USE_HFETI = true;
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unsupported FETI METHOD";
	}
	cluster.USE_KINV			= config::solver::USE_SCHUR_COMPLEMENT ? 1 : 0;
	cluster.SUBDOM_PER_CLUSTER	= number_of_subdomains_per_cluster;
	cluster.NUMBER_OF_CLUSTERS	= config::env::MPIsize;
	// ***************************************************************************************************************************

	// ***************************************************************************************************************************
	// Iter Solver Set-up
	solver.CG_max_iter	 = config::solver::ITERATIONS;
	solver.USE_GGtINV	 = 1;
	solver.epsilon		 = config::solver::EPSILON;
	solver.USE_PREC		 = config::solver::PRECONDITIONER;

	solver.USE_HFETI	 = cluster.USE_HFETI;
	solver.USE_KINV		 = cluster.USE_KINV;
	solver.USE_DYNAMIC	 = cluster.USE_DYNAMIC;
	// ***************************************************************************************************************************

	int solv_num_procs = Esutils::getEnv<int>("SOLVER_NUM_THREADS");
	int par_num_procs = Esutils::getEnv<int>("PAR_NUM_THREADS");

	cluster.PAR_NUM_THREADS	= par_num_procs;
	cluster.SOLVER_NUM_THREADS = solv_num_procs;

	solver.PAR_NUM_THREADS = par_num_procs;
	solver.SOLVER_NUM_THREADS = solv_num_procs;

	//mkl_cbwr_set(MKL_CBWR_COMPATIBLE);

}

// TODO: const parameters
void LinearSolver::init(const std::vector<int> &neighbours)
{

	// Kill Cilk+ threads
//	__cilkrts_end_cilk();

	number_of_subdomains_per_cluster = physics.K.size();

	// Overall Linear Solver Time measurement structure
	timeEvalMain.totalTime.startWithBarrier();

	TimeEvent timeSetClust(string("Solver - Set cluster")); timeSetClust.start();
	// Setup Cluster and Solver

	std::vector <eslocal> domain_list (number_of_subdomains_per_cluster, 0);
	for (eslocal i = 0; i < number_of_subdomains_per_cluster; i++) {
		domain_list[i] = i;
	}

	cluster.cluster_global_index = config::env::MPIrank + 1;
	cluster.InitClusterPC(&domain_list[0], number_of_subdomains_per_cluster);
	cluster.my_neighs = std::vector<eslocal>(neighbours.begin(), neighbours.end());
	cluster.mtype = physics.mtype;
	switch (physics.mtype) {
	case SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE:
		cluster.SYMMETRIC_SYSTEM = true;
		break;
	case SparseMatrix::MatrixType::REAL_SYMMETRIC_INDEFINITE:
		cluster.SYMMETRIC_SYSTEM = true;
		break;
	case SparseMatrix::MatrixType::REAL_UNSYMMETRIC:
		cluster.SYMMETRIC_SYSTEM = false;
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown matrix type";
	}

	if (!cluster.SYMMETRIC_SYSTEM && (config::solver::CGSOLVER != config::solver::CGSOLVERalternative::GMRES &&
      config::solver::CGSOLVER != config::solver::CGSOLVERalternative::BICGSTAB)) {
		ESINFO(GLOBAL_ERROR) << "Only GMRES or BICGSTAB solver supports non-symmetric systems.";
	}

	SINGULAR = physics.singular();


	vector<double> solver_parameters ( 10 );
	solver.Setup ( solver_parameters, cluster );

	// END - Setup Cluster and Solver
	timeSetClust.endWithBarrier(); timeEvalMain.addEvent(timeSetClust);


	// *** Setup B0 matrix *******************************************************************************************
	if (cluster.USE_HFETI == 1 ) {
		TimeEvent timeSetB0(string("Solver - Set B0")); timeSetB0.start();
		set_B0(constraints.B0);
		timeSetB0.endWithBarrier(); timeEvalMain.addEvent(timeSetB0);
	}
	// *** END - Setup B0 matrix *************************************************************************************


	// *** Setup B1 matrix *******************************************************************************************
	TimeEvent timeSetB1(string("Solver - Set B1")); timeSetB1.start();
	set_B1(constraints.B1, constraints.B1duplicity);
	timeSetB1.endWithBarrier(); timeEvalMain.addEvent(timeSetB1);
	// *** END - Setup B1 matrix *************************************************************************************


	// *** Setup R matrix ********************************************************************************************
	if (SINGULAR) {
		TimeEvent timeSetR(string("Solver - Set R"));
		timeSetR.start();

		// TODO: remove copying of R1
		#pragma omp parallel for
for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {

	//		physics.R1[d].GramSchmidtOrtho();

			cluster.domains[d].Kplus_R = physics.R1[d];
			cluster.domains[d].Kplus_R2 = physics.R2[d];
			cluster.domains[d].Kplus_Rb = physics.R1[d];
			cluster.domains[d].Kplus_Rb2 = physics.R2[d];
		}
//		cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//			cluster.domains[d].K = physics.K[d];
//			if (solver.USE_PREC == config::solver::PRECONDITIONERalternative::MAGIC) {
//				cluster.domains[d].Prec = cluster.domains[d].K;
//				cluster.domains[d].Prec.MatAddInPlace(physics.RegMat[d], 'N', -1);
//			}
//		}

		timeSetR.endWithBarrier();
		timeEvalMain.addEvent(timeSetR);
	}
	// *** END - Setup R matrix **************************************************************************************

	// *** Load InitialCondition for dirichelt
	TimeEvent timeSetInitialCondition(string("Solver - Set Dirichlet InitialCondition points"));
	timeSetInitialCondition.start();

	#pragma omp parallel for
for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster.domains[d].vec_c = constraints.B1c[d];
		cluster.domains[d].vec_lb = constraints.LB[d];
	}

	timeSetInitialCondition.endWithBarrier();
	timeEvalMain.addEvent(timeSetInitialCondition);
	// *** END - Load InitialCondition for dirichelt



	// *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *******************
	TimeEvent timeSolPrec(string("Solver - FETI Preprocessing")); timeSolPrec.start();

	#pragma omp parallel for
for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster.domains[d].lambda_map_sub = constraints.B1subdomainsMap[d];
	}
	Preprocessing( constraints.B1clustersMap );

	timeSolPrec.endWithBarrier(); timeEvalMain.addEvent(timeSolPrec);
	// *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************


	// **** Calculate Dirichlet Preconditioner ********************************
	if (config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::DIRICHLET ||
      config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET ) {
		TimeEvent timeDirPrec(string("Solver - Dirichlet Preconditioner calculation")); timeDirPrec.start();

		ESINFO(PROGRESS2) << "Calculate Dirichlet preconditioner";
        cluster.CreateDirichletPrec( physics );
       /* 
		cilk_for (size_t d = 0; d < physics.K.size(); d++) {
			SEQ_VECTOR <eslocal> perm_vec = cluster.domains[d].B1t_Dir_perm_vec;
			SEQ_VECTOR <eslocal> perm_vec_full ( physics.K[d].rows );
			SEQ_VECTOR <eslocal> perm_vec_diff ( physics.K[d].rows );

			SEQ_VECTOR <eslocal> I_row_indices_p (physics.K[d].nnz);
			SEQ_VECTOR <eslocal> J_col_indices_p (physics.K[d].nnz);

			for (size_t i = 0; i < perm_vec.size(); i++) {
				perm_vec[i] = perm_vec[i] - 1;
			}

			for (size_t i = 0; i < perm_vec_full.size(); i++) {
				perm_vec_full[i] = i;
			}

			auto it = std::set_difference( perm_vec_full.begin(), perm_vec_full.end(), perm_vec.begin(), perm_vec.end(), perm_vec_diff.begin() );
			perm_vec_diff.resize(it - perm_vec_diff.begin());

			perm_vec_full = perm_vec_diff;
			perm_vec_full.insert(perm_vec_full.end(), perm_vec.begin(), perm_vec.end());

			SparseMatrix K_modif = physics.K[d];
      SparseMatrix RegMatCRS = physics.RegMat[d];
      RegMatCRS.ConvertToCSRwithSort(0);
			K_modif.MatAddInPlace(RegMatCRS,'N',-1);
			// K_modif.RemoveLower();

			SEQ_VECTOR <SEQ_VECTOR<eslocal >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<eslocal >(2, 1));
			eslocal offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

			for (eslocal i = 0; i < K_modif.rows;i++){
				vec_I1_i2[i][0] = perm_vec_full[i];
				vec_I1_i2[i][1] = i; // position to create reverse permutation
			}

			std::sort(vec_I1_i2.begin(), vec_I1_i2.end(), [](const SEQ_VECTOR <eslocal >& a, const SEQ_VECTOR<eslocal>& b) { return a[0] < b[0]; });

			// permutations made on matrix in COO format
      K_modif.ConvertToCOO(0);
      eslocal I_index,J_index;
      bool unsymmetric=!cluster.SYMMETRIC_SYSTEM;
      for (eslocal i = 0;i<K_modif.nnz;i++){
        I_index = vec_I1_i2[K_modif.I_row_indices[i]-offset][1]+offset;
        J_index = vec_I1_i2[K_modif.J_col_indices[i]-offset][1]+offset;
        if (unsymmetric || I_index<=J_index){
          I_row_indices_p[i]=I_index;
          J_col_indices_p[i]=J_index;
        }
        else{
          I_row_indices_p[i]=J_index;
          J_col_indices_p[i]=I_index;
        }
      }
      for (eslocal i = 0; i<K_modif.nnz;i++){
        K_modif.I_row_indices[i] = I_row_indices_p[i];
        K_modif.J_col_indices[i] = J_col_indices_p[i];
      }
      K_modif.ConvertToCSRwithSort(1);
      {
			if (config::info::PRINT_MATRICES) {
				std::ofstream osS(Logging::prepareFile(d, "K_modif"));
				osS << K_modif;
				osS.close();
			}
      }


      // ------------------------------------------------------------------------------------------------------------------
      bool diagonalized_K_rr = config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET;
      //        PRECONDITIONER==NONE              - 0
      //        PRECONDITIONER==LUMPED            - 1
      //        PRECONDITIONER==WEIGHT_FUNCTION   - 2
      //        PRECONDITIONER==DIRICHLET         - 3
      //        PRECONDITIONER==SUPER_DIRICHLET   - 4
      //        
      //        When next line is uncomment, var. PRECONDITIONER==DIRICHLET and PRECONDITIONER==SUPER_DIRICHLET provide identical preconditioner.
      //        bool diagonalized_K_rr = false
      // ------------------------------------------------------------------------------------------------------------------

			eslocal sc_size = perm_vec.size();

			if (sc_size == physics.K[d].rows) {
				cluster.domains[d].Prec = physics.K[d];
				cluster.domains[d].Prec.ConvertCSRToDense(1);
        // if physics.K[d] does not contain inner DOF
			} else {

        if (config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::DIRICHLET) {
          SparseSolverCPU createSchur;
//          createSchur.msglvl=1;
          eslocal sc_size = perm_vec.size();
          createSchur.ImportMatrix_wo_Copy(K_modif);
          createSchur.Create_SC(cluster.domains[d].Prec, sc_size,false);
				  cluster.domains[d].Prec.ConvertCSRToDense(1);
        }
        else
        {
          SparseMatrix K_rr;
          SparseMatrix K_rs;
          SparseMatrix K_sr;
          SparseMatrix KsrInvKrrKrs; 

          eslocal i_start = 0;
          eslocal nonsing_size = K_modif.rows - sc_size - i_start;
          eslocal j_start = nonsing_size;

          K_rs.getSubBlockmatrix_rs(K_modif,K_rs,i_start, nonsing_size,j_start,sc_size);

          if (cluster.SYMMETRIC_SYSTEM){
            K_rs.MatTranspose(K_sr);
          }
          else
          {
            K_sr.getSubBlockmatrix_rs(K_modif,K_sr,j_start,sc_size,i_start, nonsing_size);
          }
      
          cluster.domains[d].Prec.getSubDiagBlockmatrix(K_modif,cluster.domains[d].Prec,nonsing_size,sc_size);
          SEQ_VECTOR <double> diagonals;
          SparseSolverCPU K_rr_solver;

          // K_rs is replaced by:
          // a) K_rs = 1/diag(K_rr) * K_rs          (simplified Dirichlet precond.)
          // b) K_rs =    inv(K_rr) * K_rs          (classical Dirichlet precond. assembled by own - not via PardisoSC routine)
          if (diagonalized_K_rr){
            diagonals = K_modif.getDiagonal();
            // diagonals is obtained directly from K_modif (not from K_rr to avoid assembling) thanks to its structure
            //      K_modif = [K_rr, K_rs]
            //                [K_sr, K_ss]
            // 
            for (eslocal i = 0; i < K_rs.rows; i++) {
              for (eslocal j = K_rs.CSR_I_row_indices[i]; j < K_rs.CSR_I_row_indices[i + 1]; j++) {
                K_rs.CSR_V_values[j - offset] /= diagonals[i];
              }
            }
          }
          else
          {
            K_rr.getSubDiagBlockmatrix(K_modif,K_rr,i_start, nonsing_size);
            K_rr_solver.ImportMatrix_wo_Copy(K_rr);
//            K_rr_solver.msglvl = 1;
            K_rr_solver.SolveMat_Dense(K_rs);
          }

          KsrInvKrrKrs.MatMat(K_sr,'N',K_rs);
          cluster.domains[d].Prec.MatAddInPlace(KsrInvKrrKrs,'N',-1);
//          if (!diagonalized_K_rr){
//				    cluster.domains[d].Prec.ConvertCSRToDense(1);
//          }
        }

			}

			if (config::info::PRINT_MATRICES) {
				std::ofstream osS(Logging::prepareFile(d, "S"));
				SparseMatrix SC =  cluster.domains[d].Prec;
				if (config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::DIRICHLET){
				  SC.ConvertDenseToCSR(1);
				}
				osS << SC;
				osS.close();
			}

			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
        */
		ESINFO(PROGRESS2);

		timeDirPrec.endWithBarrier(); timeEvalMain.addEvent(timeDirPrec);
	}
	// *** END - Calculate Dirichlet Preconditioner ********************************






	// *** Load Matrix K, Regularization, Schur Complement and Factorization ******************************************************************************
	TimeEvent timeSolKproc(string("Solver - K regularization and factorization"));
	timeSolKproc.start();

	ESLOG(MEMORY) << "Before K reg. and fact. process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	TimeEvent KregMem(string("Solver - K regularization mem. [MB]")); KregMem.startWithoutBarrier( GetProcessMemory_u() );

	ESINFO(PROGRESS2) << "Make K regular";
	cluster.ImportKmatrixAndRegularize( physics.K, physics.RegMat );

	KregMem.endWithoutBarrier( GetProcessMemory_u() );
	//KregMem.printLastStatMPIPerNode();

	ESLOG(MEMORY) << "After import K process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";


	// *** Computation of the Schur Complement ***************************************************************************************
	if ( cluster.USE_KINV == 1 ) {
		TimeEvent KSCMem(string("Solver - SC asm. w PARDISO-SC mem [MB]")); KSCMem.startWithoutBarrier( GetProcessMemory_u() );
		TimeEvent timeSolSC2(string("Solver - Schur Complement asm. - using PARDISO-SC")); timeSolSC2.start();

		bool USE_FLOAT = false;
		if (config::solver::SCHUR_COMPLEMENT_PREC == config::solver::SCHUR_COMPLEMENT_PRECalternative::SINGLE ) {
			USE_FLOAT = true;
		}

		cluster.Create_SC_perDomain(USE_FLOAT);

		timeSolSC2.endWithBarrier(); timeEvalMain.addEvent(timeSolSC2);
		KSCMem.endWithoutBarrier( GetProcessMemory_u() );
		//KSCMem.printLastStatMPIPerNode();

		ESLOG(MEMORY) << "After K inv. process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	} else {
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			cluster.domains[d].isOnACC = 0;
		}
	}
	// *** END - Computation of the Schur Complement **************************************************************************


	// *** K Factorization ****************************************************************************************************
	TimeEvent KFactMem(string("Solver - K factorization mem. [MB]")); KFactMem.startWithoutBarrier( GetProcessMemory_u() );
	ESINFO(PROGRESS2) << "Factorize K";

	cluster.SetupKsolvers();

	KFactMem.endWithoutBarrier( GetProcessMemory_u() );
	//KFactMem.printLastStatMPIPerNode();

	ESLOG(MEMORY) << "After K solver setup process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	timeSolKproc.endWithBarrier();
	timeEvalMain.addEvent(timeSolKproc);
	// *** END - Load Matrix K, Regularization, Schur Complement and Factorization ******************************************************************************


	// *** Setup Hybrid FETI part of the solver ********************************************************************************
	if (cluster.USE_HFETI == 1) {
		TimeEvent timeHFETIprec(string("Solver - HFETI preprocessing"));
		timeHFETIprec.start();
		cluster.SetClusterHFETI();
		timeHFETIprec.endWithBarrier();
		timeEvalMain.addEvent(timeHFETIprec);

		ESLOG(MEMORY) << "After HFETI preprocessing process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	}
	// *** END - Setup Hybrid FETI part of the solver ********************************************************************************

	if (cluster.USE_HFETI == 1 && config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::NULL_PIVOTS) {

		TimeEvent timeSolPrec2(string("Solver - FETI Preprocessing 2")); timeSolPrec2.start();

		ESLOG(MEMORY) << "Solver Preprocessing - HFETI with regularization from K matrix";
		ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent G1_perCluster_time ("Setup G1 per Cluster time - preprocessing"); G1_perCluster_time.start();
		TimeEvent G1_perCluster_mem ("Setup G1 per Cluster mem - preprocessing"); G1_perCluster_mem.startWithoutBarrier(GetProcessMemory_u());
		cluster.Create_G_perCluster   ();
		G1_perCluster_time.end(); G1_perCluster_time.printStatMPI();
		G1_perCluster_mem.endWithoutBarrier(GetProcessMemory_u()); G1_perCluster_mem.printStatMPI();

		ESLOG(MEMORY) << "Created G1 per cluster";
		ESLOG(MEMORY) << "Before HFETI create GGt process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		TimeEvent solver_Preprocessing_time ("Setup solver.Preprocessing() - pre-processing"); solver_Preprocessing_time.start();
		solver.Preprocessing ( cluster );
		solver_Preprocessing_time.end(); solver_Preprocessing_time.printStatMPI();

		cluster.Compress_G1();

		// Cleanup - of uncessary objects
		cluster._my_lamdas_map_indices.clear();
		#pragma omp parallel for
for (size_t d = 0; d < cluster.domains.size(); d++) {
			cluster.domains[d].B1.Clear();
		}

		timeSolPrec2.endWithBarrier(); timeEvalMain.addEvent(timeSolPrec2);

		ESLOG(MEMORY) << "After HFETI full preprocess process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";
	}

	// *** Final Solver Setup after K factorization **********************************************************************************
	TimeEvent timeSolAkpl(string("Solver - Set Solver After Kplus"));
	timeSolAkpl.start();
	cluster.SetClusterPC_AfterKplus();
	timeSolAkpl.endWithBarrier();
	timeEvalMain.addEvent(timeSolAkpl);
	// *** END - Final Solver Setup after K factorization ****************************************************************************
}

void LinearSolver::Solve( std::vector < std::vector < double > >  & f_vec,
		                  std::vector < std::vector < double > >  & prim_solution) {

	 TimeEvent timeSolCG(string("Solver - CG Solver runtime"));
	 timeSolCG.start();

	if (solver.USE_DYNAMIC == 0)
		solver.Solve_singular    ( cluster, f_vec, prim_solution );
	else {
		solver.Solve_non_singular( cluster, f_vec, prim_solution );
		solver.timing.totalTime.printStatMPI();
		//solver.timing.totalTime.Reset();
	}

	if ( config::mesh::AVERAGE_EDGES || config::mesh::AVERAGE_FACES ) {
		#pragma omp parallel for
for (size_t d = 0; d < cluster.domains.size(); d++) {
			vector < double >  tmp;
			tmp = prim_solution[d];
			cluster.domains[d].T.MatVec(tmp, prim_solution[d], 'N');
		}
	}

	 timeSolCG.endWithBarrier();
     timeEvalMain.addEvent(timeSolCG);

}

void LinearSolver::Postprocessing( ) {

}

void LinearSolver::finilize() {

	// Show Linear Solver Runtime Evaluation
	solver.preproc_timing.printStatsMPI();
	solver.timing.printStatsMPI();

	if (SINGULAR) solver.postproc_timing.printStatsMPI();

	solver.timeEvalAppa.printStatsMPI();

	if (SINGULAR) solver.timeEvalProj.printStatsMPI();

	if ( solver.USE_PREC != config::solver::PRECONDITIONERalternative::NONE ) solver.timeEvalPrec.printStatsMPI();

	if ( cluster.USE_HFETI == 1 ) cluster.ShowTiming();

	 timeEvalMain.totalTime.endWithBarrier();
	 timeEvalMain.printStatsMPI();
}

void LinearSolver::CheckSolution( vector < vector < double > > & prim_solution ) {
    // *** Solutin correctnes test **********************************************************************************************
	double max_v = 0.0;
		for (eslocal i = 0; i < number_of_subdomains_per_cluster; i++)
			for (size_t j = 0; j < prim_solution[i].size(); j++)
				if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );

	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.startWithoutBarrier(0.0); max_sol_ev.endWithoutBarrier(max_v);

	double max_vg;
	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	ESINFO(DETAILS) << "Maxvalue in solution = " << std::setprecision(12) << max_vg;

	//max_sol_ev.printLastStatMPIPerNode(max_vg);
	// *** END - Solutin correctnes test ******************************************************************************************

}


void LinearSolver::set_B1(
		const std::vector < SparseMatrix >			& B1_mat,
		const std::vector < std::vector <double> >	& B1_duplicity) {

	#pragma omp parallel for
for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {

		cluster.domains[d].B1 = B1_mat[d];
		cluster.domains[d].B1.type = 'G';

		cluster.domains[d].B1t = cluster.domains[d].B1;
		cluster.domains[d].B1t.MatTransposeCOO();
		cluster.domains[d].B1t.ConvertToCSRwithSort(1);

	}

	#pragma omp parallel for
for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster.domains[d].B1_scale_vec = B1_duplicity[d];
	}




}

void LinearSolver::set_B0 ( const std::vector < SparseMatrix >	& B0_mat ) {

	#pragma omp parallel for
for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster.domains[d].B0 = B0_mat[d];
		cluster.domains[d].B0.type = 'G';
		cluster.domains[d].B0.ConvertToCSRwithSort(1);
	}
}


void LinearSolver::set_R_from_K ()
{

  // Allocation of vectors on each cluster
  vector<double> norm_KR_d_pow_2; norm_KR_d_pow_2.resize(number_of_subdomains_per_cluster);
  vector<eslocal> defect_K_d; defect_K_d.resize(number_of_subdomains_per_cluster);


  // getting factors and kernels of stiffness matrix K (+ statistic)
	#pragma omp parallel for
for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {


		if (cluster.SYMMETRIC_SYSTEM) {
		  cluster.domains[d].K.get_kernel_from_K(cluster.domains[d].K,
                                            cluster.domains[d]._RegMat,
                                            cluster.domains[d].Kplus_R,
											norm_KR_d_pow_2[d],
                                            defect_K_d[d],
											d);
    }
    else
    {
		  cluster.domains[d].K.get_kernels_from_nonsym_K(cluster.domains[d].K,
                                            cluster.domains[d]._RegMat,
                                            cluster.domains[d].Kplus_R,
                                            cluster.domains[d].Kplus_R2,
											norm_KR_d_pow_2[d],
                                            defect_K_d[d],
											d);
    }


		cluster.domains[d].Kplus_Rb = cluster.domains[d].Kplus_R;
		cluster.domains[d].Kplus_Rb2 = cluster.domains[d].Kplus_R2;

	}
	ESINFO(PROGRESS2) << "K kernel detected";
  // sum of ||K*R|| (all subdomains on the cluster)
  //
  //
 // ------------------------------------ GLOBAL GETHERING OF || K * R || / || K || FIXED !!!
#ifdef VERBOSE_LEVEL_X
  double sum_per_sub_on_clst_norm_KR_d_pow_2 = 0;
	for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
     sum_per_sub_on_clst_norm_KR_d_pow_2+=norm_KR_d_pow_2[d];
  }
  //
  auto max_defect = std::max_element(defect_K_d.begin(),defect_K_d.end());
  auto min_defect = std::min_element(defect_K_d.begin(),defect_K_d.end());
  //
  auto max_norm_KR_d_pow_2 = std::max_element(norm_KR_d_pow_2.begin(),norm_KR_d_pow_2.end());
  auto min_norm_KR_d_pow_2 = std::min_element(norm_KR_d_pow_2.begin(),norm_KR_d_pow_2.end());
  //
  double MinMaxMeanNorm_MinMaxDefectNsubs[6]= {*min_norm_KR_d_pow_2, *max_norm_KR_d_pow_2,
                                           sum_per_sub_on_clst_norm_KR_d_pow_2,
                                           *min_defect,*max_defect,
                                           number_of_subdomains_per_cluster};
  //
  // gathering of statistic on MPI_rank = 0
  int recv_msg_size=6;
	MPI_Status 	  mpi_stat;
  //
  if (MPI_rank>0)
  {
    MPI_Request * mpi_send_req  = new MPI_Request [1];
    MPI_Isend(&(MinMaxMeanNorm_MinMaxDefectNsubs[0]), recv_msg_size, MPI_DOUBLE, 0, 0,
                MPI_COMM_WORLD, mpi_send_req);
  }
  else
  {
  // Allocation of vectros on MPI_rank = 0
    vector<double> norm_KR_pow_2_clusters_max; norm_KR_pow_2_clusters_max.resize(MPI_size);
    vector<double> norm_KR_pow_2_clusters_min; norm_KR_pow_2_clusters_min.resize(MPI_size);
    vector<double> norm_KR_pow_2_clusters_sum; norm_KR_pow_2_clusters_sum.resize(MPI_size);
    vector<eslocal> defect_K_clusters_max; defect_K_clusters_max.resize(MPI_size);
    vector<eslocal> defect_K_clusters_min; defect_K_clusters_min.resize(MPI_size);

    norm_KR_pow_2_clusters_min[0] = MinMaxMeanNorm_MinMaxDefectNsubs[0];
    norm_KR_pow_2_clusters_max[0] = MinMaxMeanNorm_MinMaxDefectNsubs[1];
    norm_KR_pow_2_clusters_sum[0] = MinMaxMeanNorm_MinMaxDefectNsubs[2];
    defect_K_clusters_min[0]      = MinMaxMeanNorm_MinMaxDefectNsubs[3];
    defect_K_clusters_max[0]      = MinMaxMeanNorm_MinMaxDefectNsubs[4];
    eslocal numberOfAllSubdomains     = MinMaxMeanNorm_MinMaxDefectNsubs[5];


		for (eslocal i = 1; i < MPI_size; i++) {
      //TODO: Should (may) be 'cilk_for' insead of 'for' used?
	    MPI_Recv(&(MinMaxMeanNorm_MinMaxDefectNsubs[0]), recv_msg_size, MPI_DOUBLE, i, 0,
                 MPI_COMM_WORLD, &mpi_stat);
      norm_KR_pow_2_clusters_min[i] = MinMaxMeanNorm_MinMaxDefectNsubs[0];
      norm_KR_pow_2_clusters_max[i] = MinMaxMeanNorm_MinMaxDefectNsubs[1];
      norm_KR_pow_2_clusters_sum[i] = MinMaxMeanNorm_MinMaxDefectNsubs[2];
      defect_K_clusters_min[i]      = MinMaxMeanNorm_MinMaxDefectNsubs[3];
      defect_K_clusters_max[i]      = MinMaxMeanNorm_MinMaxDefectNsubs[4];
      numberOfAllSubdomains        += MinMaxMeanNorm_MinMaxDefectNsubs[5];
    }
    //
    double sum_per_clst_norm_KR_pow_2 = 0;
	  for(eslocal c = 0; c < MPI_size; c++) {
       sum_per_clst_norm_KR_pow_2+=norm_KR_pow_2_clusters_sum[c];
    }
    // from here norm_KR is not powered by 2 anymore
    double norm_KR_clusters_mean = sqrt(sum_per_clst_norm_KR_pow_2/numberOfAllSubdomains);

    //
    auto max_defect_per_clust = std::max_element(defect_K_clusters_max.begin(),defect_K_clusters_max.end());
    auto min_defect_per_clust = std::min_element(defect_K_clusters_min.begin(),defect_K_clusters_min.end());

    auto max_norm_KR_pow_2_per_clust = std::max_element(norm_KR_pow_2_clusters_max.begin(),norm_KR_pow_2_clusters_max.end());
    auto min_norm_KR_pow_2_per_clust = std::min_element(norm_KR_pow_2_clusters_min.begin(),norm_KR_pow_2_clusters_min.end());

    double _min_norm_KR_per_clust=sqrt(*min_norm_KR_pow_2_per_clust);
    double _max_norm_KR_per_clust=sqrt(*max_norm_KR_pow_2_per_clust);



    std::ofstream os ("kernel_statistic.txt");
    os.precision(17);
    os << " *******************************************************************************************************************************\n";
    os << " ********************    K E R N E L   D E T E C T I O N    V I A    S C H U R   C O M P L E M E N T    ************************\n";
    os << " *******************************************************************************************************************************\n";
    os << " Statistics for " << numberOfAllSubdomains;
    if (MPI_size==0 &&  number_of_subdomains_per_cluster==1 ){
      os << " subdomain.\n";
     }
    else
    {
      os << " subdomains.\n";
    }
    os << " defect(K)  min:max        " << *min_defect_per_clust << " : "
                                        << *max_defect_per_clust << "\n";
    os << " ||K*R||    min:max:avg    " << _min_norm_KR_per_clust << " : "
                                        << _max_norm_KR_per_clust << " : "
                                        << norm_KR_clusters_mean << "\n";
    os.close();
  }

#endif


}


void LinearSolver::set_R (
		const Mesh &mesh
)
{
	ESINFO(GLOBAL_ERROR) << "SetR -> R is computed in assembler";

//	std::vector < std::vector < std:: vector < double > > > coordinates;
//	coordinates.resize( number_of_subdomains_per_cluster );
//
//	cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		coordinates[d].resize(mesh.coordinates().localSize(d), std::vector <double> (2, 0.0));
//		for (eslocal i = 0; i < mesh.coordinates().localSize(d); i++) {
//			coordinates[d][i][0] = mesh.coordinates().get(i, d).x;
//			coordinates[d][i][1] = mesh.coordinates().get(i, d).y;
//			coordinates[d][i][2] = mesh.coordinates().get(i, d).z;
//		}
//		cluster.domains[d].CreateKplus_R( coordinates[d] );
//
//		//TODO: *** test nesymetrickeho systemu pro GGt - smazat !!
//		if (!cluster.SYMMETRIC_SYSTEM) {
//			cluster.domains[d].Kplus_R2 = cluster.domains[d].Kplus_R;
//
////			int rows = cluster.domains[d].Kplus_R2.rows;
////			int cols = cluster.domains[d].Kplus_R2.cols;
////			for (int c = 0; c < cols; c++) {
////				int s1 = c * rows;
////				int s2 = (cols - 1 - c) * rows;
////				for (int r = 0; r < rows; r++) {
////					cluster.domains[d].Kplus_R2.dense_values[s2 + r] =
////							cluster.domains[d].Kplus_R.dense_values[s1 + r];
////				}
////			}
//
//		}
//		//***
//
//		//cluster.domains[d].Kplus_Rb = cluster.domains[d].Kplus_R;
//
//	}

}

void LinearSolver::Preprocessing( std::vector < std::vector < eslocal > > & lambda_map_sub) {

	if ( ! (cluster.USE_HFETI == 1 && config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::NULL_PIVOTS )) {
		ESLOG(MEMORY) << "Solver Preprocessing";
		ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 TimeEvent G1_perCluster_time ("Setup G1 per Cluster time - preprocessing"); G1_perCluster_time.start();
		 TimeEvent G1_perCluster_mem ("Setup G1 per Cluster mem - preprocessing"); G1_perCluster_mem.startWithoutBarrier(GetProcessMemory_u());
		cluster.Create_G_perCluster   ();
		 G1_perCluster_time.end(); G1_perCluster_time.printStatMPI();
		 G1_perCluster_mem.endWithoutBarrier(GetProcessMemory_u()); G1_perCluster_mem.printStatMPI();

		ESLOG(MEMORY) << "Created G1 per cluster";
		ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 TimeEvent solver_Preprocessing_time ("Setup solver.Preprocessing() - pre-processing"); solver_Preprocessing_time.start();
		solver.Preprocessing ( cluster );
		 solver_Preprocessing_time.end(); solver_Preprocessing_time.printStatMPI();
	}

	ESLOG(MEMORY) << "Preprocessing";
	ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	 TimeEvent cluster_SetClusterPC_time ("Setup cluster.SetClusterPC() - pre-processing"); cluster_SetClusterPC_time.start();
	cluster.SetClusterPC( lambda_map_sub );
	 cluster_SetClusterPC_time.end(); cluster_SetClusterPC_time.printStatMPI();


	ESLOG(MEMORY) << "Preprocessing end";
	ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";


}
