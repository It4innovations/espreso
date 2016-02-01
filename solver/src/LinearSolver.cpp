/*
 * LinearSolver.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#include "LinearSolver.h"

LinearSolver::LinearSolver() {

}

LinearSolver::~LinearSolver() {
	// TODO Auto-generated destructor stub
}

void LinearSolver::setup( eslocal rank, eslocal size, bool IS_SINGULAR ) {

	SINGULAR 	= IS_SINGULAR;

	if ( esconfig::solver::REGULARIZATION == 0 )
  		R_from_mesh = true	;
  	else
  		R_from_mesh = false	;

	if ( esconfig::solver::KEEP_FACTORS == 0)
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

	cluster.USE_HFETI			= esconfig::solver::FETI_METHOD;
	cluster.USE_KINV			= esconfig::solver::USE_SCHUR_COMPLEMENT;
	cluster.SUBDOM_PER_CLUSTER	= number_of_subdomains_per_cluster;
	cluster.NUMBER_OF_CLUSTERS	= MPI_size;
	cluster.DOFS_PER_NODE		= DOFS_PER_NODE;
	// ***************************************************************************************************************************

	// ***************************************************************************************************************************
	// Iter Solver Set-up
	solver.CG_max_iter	 = esconfig::solver::maxIterations;
	solver.USE_GGtINV	 = 1;
	solver.epsilon		 = esconfig::solver::epsilon;
	solver.USE_PIPECG	 = esconfig::solver::CG_SOLVER;
	solver.USE_PREC		 = esconfig::solver::PRECONDITIONER;

	solver.USE_HFETI	 = cluster.USE_HFETI;
	solver.USE_KINV		 = cluster.USE_KINV;
	solver.USE_DYNAMIC	 = cluster.USE_DYNAMIC;
	// ***************************************************************************************************************************

	/* Numbers of processors, value of SOLVER_NUM_THREADS */
	int solv_num_procs;
	char * var = getenv("SOLVER_NUM_THREADS");
    if(var != NULL)
    	sscanf( var, "%d", &solv_num_procs );
	else {
    	printf("Set environment SOLVER_NUM_THREADS to 1 - number of cores");
        exit(1);
	}

	/* Numbers of processors, value of SOLVER_NUM_THREADS */
	int par_num_procs;
	char * var2 = getenv("PAR_NUM_THREADS");
    if(var != NULL)
    	sscanf( var2, "%d", &par_num_procs );
	else {
    	printf("Set environment PAR_NUM_THREADS to 1 - number of cores");
        exit(1);
	}

	cluster.PAR_NUM_THREADS	= par_num_procs;
	cluster.SOLVER_NUM_THREADS = solv_num_procs;

	solver.PAR_NUM_THREADS = par_num_procs;
	solver.SOLVER_NUM_THREADS = solv_num_procs;

	//mkl_cbwr_set(MKL_CBWR_COMPATIBLE);

}

void LinearSolver::init(
		std::vector < SparseMatrix >	& K_mat,
		std::vector < SparseMatrix >	& B1_mat,
		std::vector < SparseMatrix >	& B0_mat,

		std::vector < std::vector <eslocal> >	& lambda_map_sub_B1,
		std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
		std::vector < std::vector <eslocal> >	& lambda_map_sub_clst,
		std::vector < std::vector <double> >	& B1_duplicity,

		std::vector < std::vector <double > >	& f_vec,
		std::vector < std::vector <double > >	& vec_c,

		std::vector < eslocal > & neigh_clusters)
{
	number_of_subdomains_per_cluster = K_mat.size();

	// Overal Linear Solver Time measurement structure
	 timeEvalMain.SetName(string("ESPRESO Linear Solver Overal Timing"));
	 timeEvalMain.totalTime.AddStartWithBarrier();


	 TimeEvent timeSetClust(string("Solver - Set cluster"));
	 timeSetClust.AddStart();

	// Setup Cluster
	std::vector <eslocal> domain_list (number_of_subdomains_per_cluster,0);
	for (eslocal i = 0; i<number_of_subdomains_per_cluster; i++)
		domain_list[i] = i;

	cluster.cluster_global_index = MPI_rank + 1;
	cluster.InitClusterPC(&domain_list[0], number_of_subdomains_per_cluster);
	cluster.my_neighs = neigh_clusters;

	// END - Setup cluster



	vector<double> solver_parameters ( 10 );
	solver.Setup ( solver_parameters, cluster );

	 timeSetClust.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetClust);


	// *** Setup B0 matrix *******************************************************************************************
	if (cluster.USE_HFETI == 1 ) {
		  TimeEvent timeSetB0(string("Solver - Set B0")); timeSetB0.AddStart();
		 set_B0(B0_mat);
		  timeSetB0.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSetB0);

	}
	// *** END - Setup B0 matrix *************************************************************************************


	// *** Setup B1 matrix *******************************************************************************************
	 TimeEvent timeSetB1(string("Solver - Set B1")); timeSetB1.AddStart();
	set_B1(B1_mat,B1_duplicity);
	 timeSetB1.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSetB1);
	// *** END - Setup B1 matrix *************************************************************************************


	// *** Setup R matrix ********************************************************************************************
	TimeEvent timeSetR(string("Solver - Set R")); timeSetR.AddStart();

	  cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		  cluster.domains[d].K = K_mat[d];
		  if ( cluster.domains[d].K.type == 'G' )
			cluster.domains[d].K.RemoveLower();
		  if ( solver.USE_PREC == 1 )
			cluster.domains[d].Prec = cluster.domains[d].K;
	  }
	  set_R_from_K();

	timeSetR.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSetR);
	// *** END - Setup R matrix **************************************************************************************


	// *** Load RHS and fix points for K regularization **************************************************************
	 TimeEvent timeSetRHS(string("Solver - Set RHS and Fix points"));
	 timeSetRHS.AddStart();

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].f = f_vec[d];

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].vec_c = vec_c[d];


	 timeSetRHS.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetRHS);
	// *** END - Load RHS and fix points for K regularization ********************************************************


	// *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *******************
	 TimeEvent timeSolPrec(string("Solver - FETI Preprocessing")); timeSolPrec.AddStart();

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].lambda_map_sub = lambda_map_sub_B1[d];

	Preprocessing( lambda_map_sub_clst );

	 timeSolPrec.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSolPrec);
	// *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************


	// *** Load Matrix K and regularization ******************************************************************************
	 TimeEvent timeSolKproc(string("Solver - K regularization and factorization"));
	 timeSolKproc.AddStart();
	if (MPI_rank == 0) std::cout << "K regularization and factorization ... " << std::endl ;
	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		if (MPI_rank == 0) std::cout << d << " " ;
		if ( d == 0 && cluster.cluster_global_index == 1) cluster.domains[d].Kplus.msglvl=1;

		if (R_from_mesh) {

			cluster.domains[d].K = K_mat[d];

			if ( cluster.domains[d].K.type == 'G' )
				cluster.domains[d].K.RemoveLower();

			if ( solver.USE_PREC == 1 )
				cluster.domains[d].Prec = cluster.domains[d].K;

			cluster.domains[d].K_regularizationFromR( cluster.domains[d].K );

		}

		// Import of Regularized matrix K into Kplus (Sparse Solver)
		cluster.domains[d].Kplus.ImportMatrix (cluster.domains[d].K);

		if (KEEP_FACTORS) {
			cluster.domains[d].Kplus.keep_factors = true;
			cluster.domains[d].Kplus.Factorization ();
		} else {
			cluster.domains[d].Kplus.keep_factors = false;
			cluster.domains[d].Kplus.MPIrank = MPI_rank;
			K_mat[d].Clear();
			cluster.domains[d].K.Clear();
		}

		cluster.domains[d].domain_prim_size = cluster.domains[d].Kplus.cols;

		if ( cluster.cluster_global_index == 1 ) { std::cout < ".";}; //{ GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }

		if ( d == 0 && cluster.cluster_global_index == 1) cluster.domains[d].Kplus.msglvl=0;
	}

	if (MPI_rank == 0) std::cout << std::endl;
	 timeSolKproc.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolKproc);
	// *** END - Load Matrix K and regularization ******************************************************************************


	// *** Setup Hybrid FETI part of the solver ********************************************************************************
	if (cluster.USE_HFETI == 1) {
		 TimeEvent timeHFETIprec(string("Solver - HFETI preprocessing"));
		 timeHFETIprec.AddStart();
		cluster.SetClusterHFETI( R_from_mesh );
		 timeHFETIprec.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeHFETIprec);
	}
	// *** END - Setup Hybrid FETI part of the solver ********************************************************************************

	// *** Computation of the Schur Complement ***************************************************************************************
	if ( cluster.USE_KINV == 1 ) {
		 TimeEvent timeSolSC1(string("Solver - Schur Complement asm. - using many RSH "));
		 timeSolSC1.AddStart();
	//	cluster.Create_Kinv_perDomain();
		 timeSolSC1.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeSolSC1);

		 TimeEvent timeSolSC2(string("Solver - Schur Complement asm. - using PARDISO-SC"));
		 timeSolSC2.AddStart();
		bool USE_FLOAT;
		cluster.Create_SC_perDomain( USE_FLOAT );
		 timeSolSC2.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeSolSC2);
	}
	// *** END - Computation of the Schur Complement *********************************************************************************


	// *** Final Solver Setup after K factorization **********************************************************************************
	 TimeEvent timeSolAkpl(string("Solver - Set Solver After Kplus"));
	 timeSolAkpl.AddStart();
	cluster.SetClusterPC_AfterKplus();
	 timeSolAkpl.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolAkpl);

}

void LinearSolver::init(
		const mesh::Mesh &mesh,

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

		std::vector < eslocal > & neigh_clusters

) {


	number_of_subdomains_per_cluster = K_mat.size();

    // Overal Linear Solver Time measurement structure
     timeEvalMain.SetName(string("ESPRESO Linear Solver Overal Timing"));
	 timeEvalMain.totalTime.AddStartWithBarrier();


	 TimeEvent timeSetClust(string("Solver - Set cluster"));
	 timeSetClust.AddStart();

	// Setup Cluster
	std::vector <eslocal> domain_list (number_of_subdomains_per_cluster,0);
	for (eslocal i = 0; i<number_of_subdomains_per_cluster; i++)
		domain_list[i] = i;

	cluster.cluster_global_index = MPI_rank + 1;
	cluster.InitClusterPC(&domain_list[0], number_of_subdomains_per_cluster);
	cluster.my_neighs = neigh_clusters;

	// END - Setup cluster



	vector<double> solver_parameters ( 10 );
	solver.Setup ( solver_parameters, cluster );

	 timeSetClust.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetClust);


	// *** Setup B0 matrix *******************************************************************************************

//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//
//		SparseMatrix B0tm;
//		B0tm = B0_mat[d];
//		B0tm.type = 'G';
//		//B0tm.PrintMatSize("B0");
//		B0tm.ConvertToCSRwithSort(1);
//		//B0tm.PrintMatSize("B0");
//
//		std::stringstream ss;
//		ss << "Bo" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		os << B0tm;
//		os.close();
//	}


	if (cluster.USE_HFETI == 1 ) {
		  TimeEvent timeSetB0(string("Solver - Set B0")); timeSetB0.AddStart();
		 set_B0(B0_mat);
		  timeSetB0.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSetB0);

	}
	// *** END - Setup B0 matrix *************************************************************************************


	// *** Setup B1 matrix *******************************************************************************************

//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//
//		SparseMatrix B0tm;
//		B0tm = B1_mat[d];
//		B0tm.type = 'G';
//		//B0tm.PrintMatSize("B0");
//		B0tm.ConvertToCSRwithSort(1);
//		//B0tm.PrintMatSize("B0");
//
//		std::stringstream ss;
//		ss << "B" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		os << B0tm;
//		os.close();
//	}

	 TimeEvent timeSetB1(string("Solver - Set B1")); timeSetB1.AddStart();
	set_B1(B1_mat,B1_duplicity);
	 timeSetB1.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSetB1);
	// *** END - Setup B1 matrix *************************************************************************************


	// *** Setup R matrix ********************************************************************************************
	TimeEvent timeSetR(string("Solver - Set R")); timeSetR.AddStart();
   if (R_from_mesh){
	    set_R(mesh);
   }
   else{
	  cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		  cluster.domains[d].K = K_mat[d];
		  if ( cluster.domains[d].K.type == 'G' )
		  	cluster.domains[d].K.RemoveLower();
		  if ( solver.USE_PREC == 1 )
		  	cluster.domains[d].Prec = cluster.domains[d].K;
	  }
      set_R_from_K();
   }


//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		std::stringstream ss;
//		ss << "R" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		SparseMatrix s = cluster.domains[d].Kplus_R;
//		s.ConvertDenseToCSR(1);
//		os << s;
//		os.close();
//	}

	timeSetR.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSetR);
	// *** END - Setup R matrix **************************************************************************************

//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		std::stringstream ss;
//		ss << "Ko" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		os << K_mat[d];
//		os.close();
//	}


	if ( esconfig::mesh::averageEdges || esconfig::mesh::averageFaces ) {
		cilk_for (int d = 0; d < T_mat.size(); d++) {

			//SpyText(T_mat[d]);

			SparseSolver Tinv;
			Tinv.mtype = 11;
			Tinv.ImportMatrix(T_mat[d]);
			Tinv.Factorization();

			//SpyText( T_mat[d] );

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

//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		std::stringstream ss;
//		ss << "KT" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		os << K_mat[d];
//		os.close();
//	}

//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		std::stringstream ss;
//		ss << "RT" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		SparseMatrix s = cluster.domains[d].Kplus_R;
//		s.ConvertDenseToCSR(1);
//		os << s;
//		os.close();
//	}

	// *** Load RHS and fix points for K regularization **************************************************************
	 TimeEvent timeSetRHS(string("Solver - Set RHS and Fix points"));
	 timeSetRHS.AddStart();

//	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
//		cluster.domains[d].f = f_vec[d];

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].vec_c = vec_c[d];




	 timeSetRHS.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSetRHS);
	// *** END - Load RHS and fix points for K regularization ********************************************************


	// *** Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *******************
	 TimeEvent timeSolPrec(string("Solver - FETI Preprocessing")); timeSolPrec.AddStart();

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].lambda_map_sub = lambda_map_sub_B1[d];

	Preprocessing( lambda_map_sub_clst );

	 timeSolPrec.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSolPrec);
	// *** END - Set up solver, create G1 per cluster, global G1, GGt, distribute GGt, factorization of GGt, compression of vector and matrices B1 and G1 *************


	// *** Load Matrix K and regularization ******************************************************************************
	 TimeEvent timeSolKproc(string("Solver - K regularization and factorization"));
	 timeSolKproc.AddStart();

	if ( cluster.cluster_global_index == 1 ) { GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }

	 TimeEvent KregMem(string("Solver - K regularization mem. [MB]")); KregMem.AddStartWOBarrier( GetProcessMemory_u() );

	if (MPI_rank == 0) std::cout << std::endl << "K regularization : ";
	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		//if (MPI_rank == 0) std::cout << "."; //<< d << " " ;
		if ( d == 0 && cluster.cluster_global_index == 1) cluster.domains[d].Kplus.msglvl=1;

	    if (R_from_mesh) {

			cluster.domains[d].K.swap(K_mat[d]);

      		if ( cluster.domains[d].K.type == 'G' )
		  		cluster.domains[d].K.RemoveLower();

//		  	if ( solver.USE_PREC == 1 )
//		  		cluster.domains[d].Prec = cluster.domains[d].K;

   			for (eslocal i = 0; i < fix_nodes[d].size(); i++)
   	 			for (eslocal d_i = 0; d_i < DOFS_PER_NODE; d_i++)
   					cluster.domains[d].fix_dofs.push_back( DOFS_PER_NODE * fix_nodes[d][i] + d_i);

			cluster.domains[d].K_regularizationFromR ( cluster.domains[d].K );

			std::vector <eslocal> ().swap (cluster.domains[d].fix_dofs);

	    }

	    if (MPI_rank == 0) std::cout << ".";

	}
	if (MPI_rank == 0) std::cout << std::endl;
	 KregMem.AddEndWOBarrier( GetProcessMemory_u() );
	 KregMem.PrintLastStatMPI_PerNode( 0.0 );

	if ( cluster.cluster_global_index == 1 ) { GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }

//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		std::stringstream ss;
//		ss << "Kreg" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		SparseMatrix s = cluster.domains[d].K;
//		os << s;
//		os.close();
//	}


//	cilk_for (int d = 0; d < T_mat.size(); d++) {
//			SparseSolver Tinv;
//			Tinv.mtype = 11;
//			Tinv.ImportMatrix(T_mat[d]);
//			Tinv.Factorization();
//
//			//SpyText( T_mat[d] );
//
//			cluster.domains[d].Kplus_R.ConvertDenseToCSR(1);
//			Tinv.SolveMat_Dense( cluster.domains[d].Kplus_R );
//			cluster.domains[d].Kplus_R.ConvertCSRToDense(1);
//	}
//
//    for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
//		std::stringstream ss;
//		ss << "RT" << d << ".txt";
//		std::ofstream os(ss.str().c_str());
//		SparseMatrix s = cluster.domains[d].Kplus_R;
//		s.ConvertDenseToCSR(1);
//		os << s;
//		os.close();
//	}


	 TimeEvent KFactMem(string("Solver - K factorization mem. [MB]")); KFactMem.AddStartWOBarrier( GetProcessMemory_u() );

	if (MPI_rank == 0) std::cout << std::endl << "K factorization : ";
	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		// Import of Regularized matrix K into Kplus (Sparse Solver)
	    cluster.domains[d].Kplus.ImportMatrix_wo_Copy (cluster.domains[d].K);

		if (KEEP_FACTORS) {
			cluster.domains[d].Kplus.keep_factors = true;
			cluster.domains[d].Kplus.Factorization ();
		} else {
			cluster.domains[d].Kplus.keep_factors = false;
			cluster.domains[d].Kplus.MPIrank = MPI_rank;
		}

		cluster.domains[d].domain_prim_size = cluster.domains[d].Kplus.cols;

		//if ( cluster.cluster_global_index == 1 ) { std::cout < ".";}; //{ GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }
	    if (MPI_rank == 0) std::cout << ".";
		if ( d == 0 && cluster.cluster_global_index == 1) cluster.domains[d].Kplus.msglvl=0;
	}
	if (MPI_rank == 0) std::cout << std::endl;
	 KFactMem.AddEndWOBarrier( GetProcessMemory_u() );
	 KFactMem.PrintLastStatMPI_PerNode( 0.0 );

	if ( cluster.cluster_global_index == 1 ) { std::cout << std::endl; GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }


	if (MPI_rank == 0) std::cout << std::endl;
	 timeSolKproc.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolKproc);
	// *** END - Load Matrix K and regularization ******************************************************************************


	// *** Setup Hybrid FETI part of the solver ********************************************************************************
	if (cluster.USE_HFETI == 1) {
		 TimeEvent timeHFETIprec(string("Solver - HFETI preprocessing"));
		 timeHFETIprec.AddStart();
		cluster.SetClusterHFETI( R_from_mesh );
		 timeHFETIprec.AddEndWithBarrier();
		 timeEvalMain.AddEvent(timeHFETIprec);

			if ( cluster.cluster_global_index == 1 ) { std::cout << std::endl; GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }

	}
	// *** END - Setup Hybrid FETI part of the solver ********************************************************************************
    //cluster.Create_G1_perCluster();

    if (cluster.USE_HFETI == 1 && !R_from_mesh ) {
    	solver.Preprocessing ( cluster );
    	cluster.G1.Clear();
    	if ( cluster.cluster_global_index == 1 ) { std::cout << std::endl; GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }

    }
    //for (int d = 0; d < cluster.domains.size(); d++)
    //	cluster.domains[d].Kplus_R = cluster.domains[d].Kplus_Rb;

	// *** Computation of the Schur Complement ***************************************************************************************
	if ( cluster.USE_KINV == 1 ) {
//		 TimeEvent timeSolSC1(string("Solver - Schur Complement asm. - using many RSH "));
//		 timeSolSC1.AddStart();
//	//	cluster.Create_Kinv_perDomain();
//		 timeSolSC1.AddEndWithBarrier();
//		 timeEvalMain.AddEvent(timeSolSC1);

		 TimeEvent KSCMem(string("Solver - SC asm. w PARDISO-SC mem [MB]")); KSCMem.AddStartWOBarrier( GetProcessMemory_u() );
		 TimeEvent timeSolSC2(string("Solver - Schur Complement asm. - using PARDISO-SC")); timeSolSC2.AddStart();
		bool USE_FLOAT = false;
		cluster.Create_SC_perDomain(USE_FLOAT);
		 timeSolSC2.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSolSC2);
		 KSCMem.AddEndWOBarrier( GetProcessMemory_u() );
		 KSCMem.PrintLastStatMPI_PerNode( 0.0 );

		if ( cluster.cluster_global_index == 1 ) { std::cout << std::endl; GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }

	}
	// *** END - Computation of the Schur Complement *********************************************************************************


	// *** Final Solver Setup after K factorization **********************************************************************************
	 TimeEvent timeSolAkpl(string("Solver - Set Solver After Kplus"));
	 timeSolAkpl.AddStart();
	cluster.SetClusterPC_AfterKplus();
	 timeSolAkpl.AddEndWithBarrier();
	 timeEvalMain.AddEvent(timeSolAkpl);
	// *** END - Final Solver Setup after K factorization ****************************************************************************





//	 timeEvalMain.totalTime.AddEndWithBarrier();
//	 timeEvalMain.PrintStatsMPI();

}

void LinearSolver::Solve( std::vector < std::vector < double > >  & f_vec,
		                  std::vector < std::vector < double > >  & prim_solution) {

	 TimeEvent timeSolCG(string("Solver - CG Solver runtime"));
	 timeSolCG.AddStart();

	if (solver.USE_DYNAMIC == 0)
		solver.Solve_singular    ( cluster, f_vec, prim_solution );
	else {
		solver.Solve_non_singular( cluster, f_vec, prim_solution );
		solver.timing.totalTime.PrintStatMPI(0.0);
		//solver.timing.totalTime.Reset();
	}

	if ( esconfig::mesh::averageEdges || esconfig::mesh::averageFaces ) {
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			vector < double >  tmp;
			tmp = prim_solution[d];
			cluster.domains[d].T.MatVec(tmp, prim_solution[d], 'N');


//			std::stringstream ss;
//			ss.precision(40);
//			ss << "sol" << d << ".txt";
//			std::ofstream os(ss.str().c_str());
//			os.precision(40);
//			os << prim_solution[d];
//			os.close();


		}
	}

	 timeSolCG.AddEndWithBarrier();
     timeEvalMain.AddEvent(timeSolCG);

}

void LinearSolver::Postprocessing( ) {

}

void LinearSolver::finilize() {

	// Show Linear Solver Runtime Evaluation
	solver.preproc_timing.PrintStatsMPI();
	solver.timing.PrintStatsMPI();

	if (SINGULAR) solver.postproc_timing.PrintStatsMPI();

	solver.timeEvalAppa.PrintStatsMPI();

	if (SINGULAR) solver.timeEvalProj.PrintStatsMPI();

	if ( solver.USE_PREC   > 0 ) solver.timeEvalPrec.PrintStatsMPI();

	if ( cluster.USE_HFETI == 1 ) cluster.ShowTiming();

	 timeEvalMain.totalTime.AddEndWithBarrier();
	 timeEvalMain.PrintStatsMPI();
}

void LinearSolver::CheckSolution( vector < vector < double > > & prim_solution ) {
    // *** Solutin correctnes test **********************************************************************************************
	double max_v = 0.0;
		for (eslocal i = 0; i < number_of_subdomains_per_cluster; i++)
			for (eslocal j = 0; j < prim_solution[i].size(); j++)
				if ( fabs ( prim_solution[i][j] ) > max_v) max_v = fabs( prim_solution[i][j] );

	TimeEvent max_sol_ev ("Max solution value "); max_sol_ev.AddStartWOBarrier(0.0); max_sol_ev.AddEndWOBarrier(max_v);

	std::cout.precision(12);
	double max_vg;
	MPI_Reduce(&max_v, &max_vg, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
	if (MPI_rank == 0)
		std::cout << " Max value in_solution = " << max_vg << std::endl;

	max_sol_ev.PrintLastStatMPI_PerNode(max_vg);
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

  // Allocation of vectros on each cluster
  vector<double> norm_KR_d_pow_2; norm_KR_d_pow_2.resize(number_of_subdomains_per_cluster);
  vector<eslocal> defect_K_d; defect_K_d.resize(number_of_subdomains_per_cluster);

  
  // getting factors and kernels of stiffness matrix K (+ statistic)
	cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		cluster.domains[d].K.get_kernel_from_K(cluster.domains[d].K,
                                            cluster.domains[d].Kplus_R,&(norm_KR_d_pow_2[d]), 
                                            &(defect_K_d[d]));
    
		cluster.domains[d].Kplus_Rb = cluster.domains[d].Kplus_R;
	}
  // sum of ||K*R|| (all subdomains on the cluster)
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

  std::cout<<"###############################################################################################################\n";
  std::cout<<"##############    K E R N E L   D E T E C T I O N    V I A    S C H U R   C O M P L E M E N T    ##############\n";
  std::cout<<"###############################################################################################################\n";
  std::cout<< " Statistics for " << numberOfAllSubdomains  << " subdomains.\n";
  std::cout<< " defect_K    min:max        " << *min_defect_per_clust << " : "
                                             << *max_defect_per_clust << "\n"; 
  std::cout<< " norm_KR     min:max:mean   " << _min_norm_KR_per_clust << " : "
                                             << _max_norm_KR_per_clust << " : "
                                             << norm_KR_clusters_mean << "\n";
  std::cout<<"##########################################################################################################\n";

  }





}


void LinearSolver::set_R (
		const mesh::Mesh &mesh
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

	if (MPI_rank == 0) {
		cout << " ******************************************************************************************************************************* " << endl;
		cout << " *** SetSolverPreprocessing **************************************************************************************************** " << endl;

		GetProcessMemoryStat_u ( );
		GetMemoryStat_u( );
		cout << " Solver - Creating G1 per cluster ... " << endl;
	}

	TimeEvent G1_perCluster_time ("Setup G1 per Cluster time - preprocessing");
	G1_perCluster_time.AddStart(omp_get_wtime());

	TimeEvent G1_perCluster_mem ("Setup G1 per Cluster mem - preprocessing");
	G1_perCluster_mem.AddStartWOBarrier(GetProcessMemory_u());

	cluster.Create_G1_perCluster   ();

	G1_perCluster_time.AddEnd(omp_get_wtime());
	G1_perCluster_time.PrintStatMPI(0.0);

	G1_perCluster_mem.AddEndWOBarrier(GetProcessMemory_u());
	G1_perCluster_mem.PrintStatMPI(0.0);

	if (MPI_rank == 0) {
		GetProcessMemoryStat_u ( );
		GetMemoryStat_u ( );
		cout << " Solver - CreateVec_d_perCluster ... " << endl;
	}

//	TimeEvent Vec_d_perCluster_time ("Setup Vec d per Cluster - preprocessing");
//	Vec_d_perCluster_time.AddStart(omp_get_wtime());
//
//	cluster.CreateVec_d_perCluster ();
//
//	Vec_d_perCluster_time.AddEnd(omp_get_wtime());
//	Vec_d_perCluster_time.PrintStatMPI(0.0);

	if (MPI_rank == 0) {
		GetProcessMemoryStat_u ( );
		GetMemoryStat_u( );
		cout << " Solver - Creating Global G1 and Running preprocessing (create GGt) ... " << endl;
	}

	TimeEvent solver_Preprocessing_time ("Setup solver.Preprocessing() - preprocessing");
	solver_Preprocessing_time.AddStart(omp_get_wtime());

	//cluster.my_neighs = neigh_domains;

	solver.Preprocessing ( cluster );

	solver_Preprocessing_time.AddEnd(omp_get_wtime());
	solver_Preprocessing_time.PrintStatMPI(0.0);

	if (MPI_rank == 0) {
		GetProcessMemoryStat_u ( );
		GetMemoryStat_u ( );
		cout << " Solver - SetClusterPC - lambda map sub a neigh domains ... " << endl;
	}

	TimeEvent cluster_SetClusterPC_time ("Setup cluster.SetClusterPC() - preprocessing");
	cluster_SetClusterPC_time.AddStart(omp_get_wtime());

	cluster.SetClusterPC( lambda_map_sub ); // USE_DYNAMIC, USE_KINV

	cluster_SetClusterPC_time.AddEnd(omp_get_wtime());
	cluster_SetClusterPC_time.PrintStatMPI(0.0);

	if (MPI_rank == 0) {
		GetProcessMemoryStat_u ( );
		GetMemoryStat_u ( );

		cout << " *** END - SetSolverPreprocessing ********************************************************************************************** " << endl;
		cout << " ******************************************************************************************************************************* " << endl;
		cout << endl;
	}


}
