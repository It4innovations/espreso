/*
 * LinearSolver.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: lriha
 */

#include "LinearSolver.h"

void Set_CSR_Matrix   (
		SparseMatrix & Mat,
		eslocal n_rows,
		eslocal n_cols,
		eslocal * rows,
		eslocal * cols,
		double * vals,
		char type ) {

	int nnz = rows[n_rows];
	int offset = (rows[0]) ? 0 : 1;
	nnz -= rows[0];

	Mat.CSR_I_row_indices.resize(n_rows+1);
	Mat.CSR_J_col_indices.resize(nnz);
	Mat.CSR_V_values	 .resize(nnz);

	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (int i = 0; i < Mat.CSR_I_row_indices.size(); i++)
		Mat.CSR_I_row_indices[i] = rows[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (int i = 0; i < Mat.CSR_J_col_indices.size(); i++)
		Mat.CSR_J_col_indices[i] = cols[i] + offset;

	copy(vals, vals + nnz, Mat.CSR_V_values.begin());

	Mat.cols = n_cols;
	Mat.rows = n_rows;
	Mat.nnz  = nnz;
	Mat.type = type;
}

void Set_COO_Matrix   (
		SparseMatrix & Mat,
		eslocal n_rows,
		eslocal n_cols,
		eslocal nnz,
		eslocal * I_rows,
		eslocal	 * J_cols,
		double * V_vals,
		char type,
		int indexing ) {

	Mat.I_row_indices.resize(nnz);
	Mat.J_col_indices.resize(nnz);
	Mat.V_values	 .resize(nnz);
	int offset = indexing ? 0 : 1;

	//copy(rows, rows + n_cols + 1, K.CSR_I_row_indices.begin());
	for (int i = 0; i < Mat.I_row_indices.size(); i++)
		Mat.I_row_indices[i] = I_rows[i] + offset;

	//copy(cols, cols + nnz, K.CSR_J_col_indices.begin());
	for (int i = 0; i < Mat.J_col_indices.size(); i++)
		Mat.J_col_indices[i] = J_cols[i] + offset;

	copy(V_vals, V_vals + nnz, Mat.V_values.begin());

	Mat.cols = n_cols;
	Mat.rows = n_rows;
	Mat.nnz  = nnz;
	Mat.type = type;

}



LinearSolver::LinearSolver() {

}

LinearSolver::~LinearSolver() {
	// TODO Auto-generated destructor stub
}

void LinearSolver::setup( int rank, int size, bool IS_SINGULAR ) {

	SINGULAR = IS_SINGULAR;

	DOFS_PER_NODE = 3;

	KEEP_FACTORS = false; // only suported by MKL Pardiso so far

    MPI_rank = rank;
    MPI_size = size;

    // ***************************************************************************************************************************
	// Cluster structure  setup
	if ( SINGULAR )
		cluster.USE_DYNAMIC		= 0;
	else
		cluster.USE_DYNAMIC		= 1;

	cluster.USE_HFETI			= 1;
	cluster.USE_KINV			= 0;
	cluster.SUBDOM_PER_CLUSTER	= number_of_subdomains_per_cluster;
	cluster.NUMBER_OF_CLUSTERS	= MPI_size;
	cluster.DOFS_PER_NODE		= DOFS_PER_NODE;
	// ***************************************************************************************************************************

	// ***************************************************************************************************************************
	// Iter Solver Set-up
	solver.CG_max_iter	 = 100;
	solver.USE_GGtINV	 = 1;
	solver.epsilon		 = 0.0001;
	solver.USE_PIPECG	 = 0;
	solver.USE_PREC		 = 1;

	solver.USE_HFETI	 = cluster.USE_HFETI;
	solver.USE_KINV		 = cluster.USE_KINV;
	solver.USE_DYNAMIC	 = cluster.USE_DYNAMIC;
	// ***************************************************************************************************************************

}

void LinearSolver::init(
		const mesh::Mesh &mesh,

		std::vector < SparseMatrix >	& K_mat,
		std::vector < SparseMatrix >	& B1_mat,
		std::vector < SparseMatrix >	& B0_mat,

		std::vector < std::vector <eslocal> >	& lambda_map_sub_B1,
		std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
		std::vector < std::vector <eslocal> >	& lambda_map_sub_clst,
		std::vector < std::vector <double> >	& B1_duplicity,

		std::vector < std::vector <double > >	& f_vec,
		std::vector < std::vector <eslocal > >	& fix_nodes,
		std::vector < std::vector <eslocal> >	& l2g_vec,

		std::vector < eslocal > & neigh_clusters

) {

	number_of_subdomains_per_cluster = K_mat.size();


    // Overal Linear Solver Time measurement structure
     timeEvalMain.SetName(string("ESPRESO Linear Solver Overal Timing"));
	 timeEvalMain.totalTime.AddStartWithBarrier();


	 TimeEvent timeSetClust(string("Solver - Set cluster"));
	 timeSetClust.AddStart();

	// Setup Cluster
	std::vector <int> domain_list (number_of_subdomains_per_cluster,0);
	for (int i = 0; i<number_of_subdomains_per_cluster; i++)
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
	set_R(l2g_vec, mesh);
	 timeSetR.AddEndWithBarrier(); timeEvalMain.AddEvent(timeSetR);
	// *** END - Setup R matrix **************************************************************************************


	// *** Load RHS and fix points for K regularization **************************************************************
	 TimeEvent timeSetRHS(string("Solver - Set RHS and Fix points"));
	 timeSetRHS.AddStart();

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		cluster.domains[d].f = f_vec[d];

	cilk_for (eslocal d = 0; d < number_of_subdomains_per_cluster; d++)
		for (int i = 0; i < fix_nodes[d].size(); i++)
 			for (int d_i = 0; d_i < 3; d_i++)
				cluster.domains[d].fix_dofs.push_back( 3 * fix_nodes[d][i] + d_i);

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

		cluster.domains[d].K = K_mat[d];

		if ( cluster.domains[d].K.type == 'G' )
			cluster.domains[d].K.RemoveLower();

		//TODO: POZOR - zbytecne kopiruju - pokud se nepouziva LUMPED
		cluster.domains[d].Prec = cluster.domains[d].K;

		cluster.domains[d].K_regularizationFromR( );

		if (KEEP_FACTORS) {
			cluster.domains[d].Kplus.keep_factors = true;
			cluster.domains[d].Kplus.Factorization ();
		} else {
			cluster.domains[d].Kplus.keep_factors = false;
		}

		cluster.domains[d].domain_prim_size = cluster.domains[d].Kplus.cols;

		if ( cluster.cluster_global_index == 1 ) { GetMemoryStat_u ( ); GetProcessMemoryStat_u ( ); }

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
		cluster.SetClusterHFETI();
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
		cluster.Create_SC_perDomain();
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
	// *** END - Final Solver Setup after K factorization ****************************************************************************





//	 timeEvalMain.totalTime.AddEndWithBarrier();
//	 timeEvalMain.PrintStatsMPI();

}

void LinearSolver::Solve(
		std::vector < std::vector <double > >	& f_vec,
		vector < vector < double > > & prim_solution) {

	 TimeEvent timeSolCG(string("Solver - CG Solver runtime"));
	 timeSolCG.AddStart();

	if (solver.USE_DYNAMIC == 0)
		solver.Solve_singular    ( cluster, f_vec, prim_solution );
	else {
		solver.Solve_non_singular( cluster, f_vec, prim_solution );
		solver.timing.totalTime.PrintStatMPI(0.0);
		//solver.timing.totalTime.Reset();
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

	if ( solver.USE_PREC   == 1 ) solver.timeEvalPrec.PrintStatsMPI();

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
		std::vector < SparseMatrix >	& B1_mat,
		std::vector < std::vector <double> >	    & B1_duplicity) {

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

void LinearSolver::set_R (
		std::vector < std::vector <eslocal> >	& l2g_vec,
		const mesh::Mesh &mesh
)
{

	cilk_for(eslocal d = 0; d < number_of_subdomains_per_cluster; d++) {
		for (int i = 0; i < l2g_vec[d].size(); i++) {
			std::vector <double> tmp_vec (3,0);
			tmp_vec[0] = mesh.coordinates()[l2g_vec[d][i]].x;
			tmp_vec[1] = mesh.coordinates()[l2g_vec[d][i]].y;
			tmp_vec[2] = mesh.coordinates()[l2g_vec[d][i]].z;
			cluster.domains[d].coordinates.push_back(tmp_vec);
		}
		cluster.domains[d].CreateKplus_R();
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
