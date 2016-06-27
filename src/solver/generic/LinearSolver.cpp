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

LinearSolver::LinearSolver(): timeEvalMain("ESPRESO Solver Overal Timing") {

}

LinearSolver::~LinearSolver() {
	// TODO Auto-generated destructor stub
}

void LinearSolver::setup( eslocal rank, eslocal size, bool IS_SINGULAR ) {

	SINGULAR 	= IS_SINGULAR;
	R_from_mesh = config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::FIX_POINTS;

	if (!config::solver::KEEP_FACTORS)
		KEEP_FACTORS = false; // only suported by MKL Pardiso so far
	else
		KEEP_FACTORS = true;

    MPI_rank = rank;
    MPI_size = size;

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
	cluster.NUMBER_OF_CLUSTERS	= MPI_size;
	cluster.DOFS_PER_NODE		= DOFS_PER_NODE;
	// ***************************************************************************************************************************

	// ***************************************************************************************************************************
	// Iter Solver Set-up
	solver.CG_max_iter	 = config::solver::ITERATIONS;
	solver.USE_GGtINV	 = 1;
	solver.epsilon		 = config::solver::EPSILON;
	solver.USE_PIPECG	 = config::solver::CG_SOLVER == config::solver::CG_SOLVERalternative::PIPELINED;
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

void LinearSolver::init(
		const Mesh &mesh,

		std::vector < SparseMatrix >	& K_mat,
		std::vector < SparseMatrix >	& T_mat,
		std::vector < SparseMatrix >	& B1_mat,
		std::vector < SparseMatrix >	& B0_mat,

		std::vector < std::vector <eslocal> >	& lambda_map_sub_B1,
		std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
		std::vector < std::vector <eslocal> >	& lambda_map_sub_clst,
		std::vector < std::vector <double> >	& B1_duplicity,

		std::vector < std::vector <double > >	& f_vec,
		std::vector < std::vector <double > >	& vec_c,

		const std::vector < std::vector <eslocal > >	& fix_nodes,

		const std::vector < int > & neigh_clusters

) {

	number_of_subdomains_per_cluster = K_mat.size();

	// Overall Linear Solver Time measurement structure
	 timeEvalMain.totalTime.startWithBarrier();

	 TimeEvent timeSetClust(string("Solver - Set cluster")); timeSetClust.start();
	// Setup Cluster and Solver
	std::vector <eslocal> domain_list (number_of_subdomains_per_cluster,0);
	for (eslocal i = 0; i<number_of_subdomains_per_cluster; i++)
		domain_list[i] = i;

	cluster.cluster_global_index = MPI_rank + 1;
	cluster.InitClusterPC(&domain_list[0], number_of_subdomains_per_cluster);
	cluster.my_neighs = std::vector<eslocal>(neigh_clusters.begin(), neigh_clusters.end());

	vector<double> solver_parameters ( 10 );
	solver.Setup ( solver_parameters, cluster );
	// END - Setup Cluster and Solver
	 timeSetClust.endWithBarrier(); timeEvalMain.addEvent(timeSetClust);


	// *** Setup B0 matrix *******************************************************************************************
	if (cluster.USE_HFETI == 1 ) {
		 TimeEvent timeSetB0(string("Solver - Set B0")); timeSetB0.start();
		set_B0(B0_mat);
		 timeSetB0.endWithBarrier(); timeEvalMain.addEvent(timeSetB0);

	}
	// *** END - Setup B0 matrix *************************************************************************************


	// *** Setup B1 matrix *******************************************************************************************
	 TimeEvent timeSetB1(string("Solver - Set B1")); timeSetB1.start();
	set_B1(B1_mat,B1_duplicity);
	 timeSetB1.endWithBarrier(); timeEvalMain.addEvent(timeSetB1);
	// *** END - Setup B1 matrix *************************************************************************************


	// *** Setup R matrix ********************************************************************************************
	 if (SINGULAR) {
		 TimeEvent timeSetR(string("Solver - Set R"));
		 timeSetR.start();
		 if (R_from_mesh) {
			 set_R(mesh);
		 } else {
			 cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
				 cluster.domains[d].K = K_mat[d];
				 if (cluster.domains[d].K.type == 'G')
					 cluster.domains[d].K.RemoveLower();
				 if (solver.USE_PREC == config::solver::PRECONDITIONERalternative::MAGIC)
					 cluster.domains[d].Prec = cluster.domains[d].K;
			 }
			 set_R_from_K();
		 }

		 if (config::info::PRINT_MATRICES) {
			 for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
				 SparseMatrix s = cluster.domains[d].Kplus_R;
				 s.ConvertDenseToCSR(1);

				 std::ofstream os(Logging::prepareFile(d, "R"));
				 os << s;
				 os.close();
			 }
		 }

		 timeSetR.endWithBarrier();
		 timeEvalMain.addEvent(timeSetR);
		 // *** END - Setup R matrix **************************************************************************************
	 }


	// *** HTFETI - averaging objects
	if ( config::mesh::AVERAGE_EDGES || config::mesh::AVERAGE_FACES ) {
		cilk_for (int d = 0; d < T_mat.size(); d++) {

			SparseSolverCPU Tinv;
			Tinv.mtype = 11;
			Tinv.ImportMatrix(T_mat[d]);
			std::stringstream ss;
			ss << "Init averaging -> rank: " << config::env::MPIrank << ", subdomain: " << d;
			Tinv.Factorization(ss.str());

			cluster.domains[d].Kplus_R.ConvertDenseToCSR(1);
			Tinv.SolveMat_Dense( cluster.domains[d].Kplus_R );
			cluster.domains[d].Kplus_R.ConvertCSRToDense(1);

			cluster.domains[d].T = T_mat[d];

			SparseMatrix Ktmp;
			Ktmp = K_mat[d];

			Ktmp.MatMat( T_mat[d], 'T', K_mat[d] );

			K_mat[d].MatMat( Ktmp, 'N', T_mat[d] );

			vector <double > tmp;
			tmp = f_vec[d];

			T_mat[d].MatVec(tmp,f_vec[d],'T');


		}

	}
	// *** END - HTFETI - averaging objects


	if (config::info::PRINT_MATRICES) {
		for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
			SparseMatrix RT = cluster.domains[d].Kplus_R;
			RT.ConvertDenseToCSR(1);

			std::ofstream osKT(Logging::prepareFile(d, "KT"));
			osKT << K_mat[d];
			osKT.close();

			std::ofstream osRT(Logging::prepareFile(d, "RT"));
			osRT << RT;
			osRT.close();
		}
	}


	// *** Load RHS for dirichelt
	 TimeEvent timeSetRHS(string("Solver - Set Dirichlet RHS points"));
	 timeSetRHS.start();

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].vec_c = vec_c[d];

	 timeSetRHS.endWithBarrier();
	 timeEvalMain.addEvent(timeSetRHS);
	// *** END - Load RHS for dirichelt



	// *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *******************
	 TimeEvent timeSolPrec(string("Solver - FETI Preprocessing")); timeSolPrec.start();

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].lambda_map_sub = lambda_map_sub_B1[d];
	Preprocessing( lambda_map_sub_clst );

	timeSolPrec.endWithBarrier(); timeEvalMain.addEvent(timeSolPrec);
	 // *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************


	// **** Calculate Dirichlet Preconditioner ********************************
	if (config::solver::PRECONDITIONER == config::solver::PRECONDITIONERalternative::DIRICHLET ) {
		TimeEvent timeDirPrec(string("Solver - Dirichlet Preconditioner calculation")); timeDirPrec.start();

		ESINFO(PROGRESS2) << "Calculate Dirichlet preconditioner";
		cilk_for (int d = 0; d < K_mat.size(); d++) {
			SEQ_VECTOR <eslocal> perm_vec = cluster.domains[d].B1t_Dir_perm_vec;
			SEQ_VECTOR <eslocal> perm_vec_full ( K_mat[d].rows );
			SEQ_VECTOR <eslocal> perm_vec_diff ( K_mat[d].rows );

			SEQ_VECTOR <eslocal> I_row_indices_p (K_mat[d].nnz);
			SEQ_VECTOR <eslocal> J_col_indices_p (K_mat[d].nnz);

			for (eslocal i = 0; i < perm_vec.size(); i++)
				perm_vec[i] = perm_vec[i] - 1;

			for (eslocal i = 0; i < perm_vec_full.size(); i++)
				perm_vec_full[i] = i;

			auto it = std::set_difference( perm_vec_full.begin(), perm_vec_full.end(), perm_vec.begin(), perm_vec.end(), perm_vec_diff.begin() );
			perm_vec_diff.resize(it - perm_vec_diff.begin());

			perm_vec_full = perm_vec_diff;
			perm_vec_full.insert(perm_vec_full.end(), perm_vec.begin(), perm_vec.end());

			SparseMatrix K_modif = K_mat[d];
			K_modif.RemoveLower();

			SEQ_VECTOR <SEQ_VECTOR<eslocal >> vec_I1_i2(K_modif.rows, SEQ_VECTOR<eslocal >(2, 1));
			eslocal offset = K_modif.CSR_I_row_indices[0] ? 1 : 0;

			for (eslocal i = 0; i < K_modif.rows;i++){
				vec_I1_i2[i][0] = perm_vec_full[i];
				vec_I1_i2[i][1] = i; // position to create reverse permutation
			}

			std::sort(vec_I1_i2.begin(), vec_I1_i2.end(),
					[](const SEQ_VECTOR <eslocal >& a, const SEQ_VECTOR <eslocal >& b) {
				return a[0] < b[0];
			});


			// permutations made on matrix in COO format
			K_modif.ConvertToCOO(0);
			eslocal I_index,J_index;
			for (eslocal i = 0;i<K_modif.nnz;i++){
				I_index = vec_I1_i2[K_modif.I_row_indices[i]-offset][1]+offset;
				J_index = vec_I1_i2[K_modif.J_col_indices[i]-offset][1]+offset;
				if (I_index>J_index){
					I_row_indices_p[i]=J_index; J_col_indices_p[i]=I_index;
				}
				else{
					I_row_indices_p[i]=I_index; J_col_indices_p[i]=J_index;
				}
			}
			//
			for (eslocal i = 0; i<K_modif.nnz;i++){
				K_modif.I_row_indices[i] = I_row_indices_p[i];
				K_modif.J_col_indices[i] = J_col_indices_p[i];
			}
			K_modif.ConvertToCSRwithSort(1);


			eslocal SC_SIZE = perm_vec.size();
			SparseMatrix S;

			if (SC_SIZE == K_mat[d].rows){
				S = K_mat[d];
				cluster.domains[d].Prec = S;
			} else {
				SparseSolverCPU createSchur;
				eslocal SC_SIZE = perm_vec.size();
				createSchur.ImportMatrix_wo_Copy(K_modif);
				createSchur.Create_SC(S,SC_SIZE,false);
				S.type='S';

				cluster.domains[d].Prec = S;
			}
	    if (config::info::PRINT_MATRICES) {
        std::ofstream osS(Logging::prepareFile(d, "S"));
        osS << S;
        osS.close();
      }



			ESINFO(PROGRESS2) << Info::plain() << ".";
		}
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
	cluster.ImportKmatrixAndRegularize( K_mat, fix_nodes );

	 KregMem.endWithoutBarrier( GetProcessMemory_u() );
	 //KregMem.printLastStatMPIPerNode();

	ESLOG(MEMORY) << "After import K process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
	ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

	if (config::info::PRINT_MATRICES) {
		for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
			std::ofstream os(Logging::prepareFile(d, "Kreg"));
			SparseMatrix s = cluster.domains[d].K;
			os << s;
			os.close();

//			std::ofstream osK(Logging::prepareFile(d, "K").c_str());
//			osK << _K[d];
//			osK.close();
//			std::ofstream osT(Logging::prepareFile(d, "T").c_str());
//			osT << _T[d];
//			osT.close();

//		cluster.domains[d].vec_c = vec_c[d];


      std::ofstream os_vec_c(Logging::prepareFile(d, "c").c_str());
			os_vec_c << vec_c[d];
			os_vec_c.close();



      std::ofstream os_in_weight(Logging::prepareFile(d, "loc_ind_weight").c_str());
			os_in_weight << lambda_map_sub_B1[d];
			os_in_weight.close();
//      lambda_map_sub_B1[d]

			std::ofstream os_weigth(Logging::prepareFile(d, "weight").c_str());
			os_weigth << B1_duplicity[d];
			os_weigth.close();

		}
	}


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
		for (int d = 0; d < cluster.domains.size(); d++) {
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
		cluster.SetClusterHFETI( R_from_mesh );
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
		cluster.Create_G1_perCluster   ();
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
    	cilk_for (eslocal d = 0; d < cluster.domains.size(); d++)
    		cluster.domains[d].B1.Clear();

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
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
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
			for (eslocal j = 0; j < prim_solution[i].size(); j++)
				if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );

	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.startWithoutBarrier(0.0); max_sol_ev.endWithoutBarrier(max_v);

	double max_vg;
	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	ESINFO(DETAILS) << "Maxvalue in solution = " << std::setprecision(12) << max_vg;

	//max_sol_ev.printLastStatMPIPerNode(max_vg);
	// *** END - Solutin correctnes test ******************************************************************************************

}


void LinearSolver::set_B1(
		std::vector < SparseMatrix >			& B1_mat,
		std::vector < std::vector <double> >	& B1_duplicity) {

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {

		cluster.domains[d].B1 = B1_mat[d];
		cluster.domains[d].B1.type = 'G';

		cluster.domains[d].B1t = cluster.domains[d].B1;
		cluster.domains[d].B1t.MatTransposeCOO();
		cluster.domains[d].B1t.ConvertToCSRwithSort(1);

	}

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster.domains[d].B1_scale_vec = B1_duplicity[d];
	}




}

void LinearSolver::set_B0 ( std::vector < SparseMatrix >	& B0_mat ) {

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
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
	cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster.domains[d].K.get_kernel_from_K(cluster.domains[d].K,
                                            cluster.domains[d]._RegMat,
                                            cluster.domains[d].Kplus_R,&(norm_KR_d_pow_2[d]),
                                            &(defect_K_d[d]),d);


		cluster.domains[d].Kplus_Rb = cluster.domains[d].Kplus_R;
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
	std::vector < std::vector < std:: vector < double > > > coordinates;
	coordinates.resize( number_of_subdomains_per_cluster );

	cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		coordinates[d].resize(mesh.coordinates().localSize(d), std::vector <double> (3, 0.0));
		for (eslocal i = 0; i < mesh.coordinates().localSize(d); i++) {
			coordinates[d][i][0] = mesh.coordinates().get(i, d).x;
			coordinates[d][i][1] = mesh.coordinates().get(i, d).y;
			coordinates[d][i][2] = mesh.coordinates().get(i, d).z;
		}
		cluster.domains[d].CreateKplus_R( coordinates[d] );
		//cluster.domains[d].Kplus_Rb = cluster.domains[d].Kplus_R;

	}

}

void LinearSolver::Preprocessing( std::vector < std::vector < eslocal > > & lambda_map_sub) {

	if ( ! (cluster.USE_HFETI == 1 && config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::NULL_PIVOTS )) {
		ESLOG(MEMORY) << "Solver Preprocessing";
		ESLOG(MEMORY) << "process " << config::env::MPIrank << " uses " << Measure::processMemory() << " MB";
		ESLOG(MEMORY) << "Total used RAM " << Measure::usedRAM() << "/" << Measure::availableRAM() << " [MB]";

		 TimeEvent G1_perCluster_time ("Setup G1 per Cluster time - preprocessing"); G1_perCluster_time.start();
		 TimeEvent G1_perCluster_mem ("Setup G1 per Cluster mem - preprocessing"); G1_perCluster_mem.startWithoutBarrier(GetProcessMemory_u());
		cluster.Create_G1_perCluster   ();
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
