//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#include "IterSolver.h"

IterSolver::IterSolver()
{

	// Timing objects

	// Main timing object for main CG loop
	timing.			SetName("Main CG loop timing ");
	preproc_timing.	SetName("Preprocessing timing ");


	timeEvalAppa.	SetName(string("Apply Kplus timing "));
	apa_B1t.		SetName(string("x = B1t * lambda "));
	apa_kplus.		SetName(string("multKplus(local or global) "));
	apa_B1.			SetName(string("lambda = B1 * x "));
	apa_allred.		SetName(string("All_Reduce_lambdas "));
	timeEvalAppa.AddEvent(apa_B1t);
	timeEvalAppa.AddEvent(apa_kplus);
	timeEvalAppa.AddEvent(apa_B1);
	timeEvalAppa.AddEvent(apa_allred);

	timeEvalPrec.	SetName(string("Apply Precond. timing "));
	prec_kplus.		SetName(string("B1 * P * B1t "));
	prec_allred.	SetName(string("All_Reduce_lambdas "));
	timeEvalPrec.AddEvent(prec_kplus);
	timeEvalPrec.AddEvent(prec_allred);


	timeEvalProj.	SetName(string("Projector timing "));
	proj_G1t.		SetName(string("x = G1 * lambda "));
	proj_Gthr.		SetName(string("MPI_gather - collective "));
	proj_GGt.		SetName(string("GGt Solve on master node "));
	proj_Sctr.		SetName(string("MPI_Scatter - collective "));
	proj_Gx.		SetName(string("lambda = G1t * x "));
	proj_allred.	SetName(string("All_Reduce_lambdas "));
	timeEvalProj.AddEvent(proj_G1t);
	timeEvalProj.AddEvent(proj_Gthr);
	timeEvalProj.AddEvent(proj_GGt);
	timeEvalProj.AddEvent(proj_Sctr);
	timeEvalProj.AddEvent(proj_Gx);
	timeEvalProj.AddEvent(proj_allred);

	ddot_time.	SetName(string("Parallel DDOT - alpha and gamma"));
	proj_time.	SetName(string("Projector_l "));
	appA_time.	SetName(string("ApplyA_l "));
	vec_time.	SetName(string("vector processing in CG "));
	norm_time.	SetName(string("parallel DDOT - norm "));

	proj1_time.	SetName(string("Projector_l - before PREC "));
	proj2_time.	SetName(string("Projector_l - after PREC "));
	prec_time.	SetName(string("Preconditioner "));
	ddot_alpha.	SetName(string("2x ddot for Alpha "));
	ddot_beta.	SetName(string("2x ddot for Beta "));



}


IterSolver::~IterSolver()
{

}

void IterSolver::Setup ( SEQ_VECTOR <double> & parameters , Cluster & cluster_in )
{

	// *** MPI variables  **********************************************************
	//mpi_rank;	//mpi_root;	//mpi_size;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	/* get number of processes */
	mpi_root = 0;

	//CG_max_iter            = 5000;
	NumberOfTimeIterations = 100;

}

void IterSolver::Preprocessing ( Cluster & cluster )
{

	preproc_timing.totalTime.AddStart(omp_get_wtime());

	// ****************************************************************************
	// *** Coarse problem - Make GGt **********************************************

	TimeEvent createGGT_time("Time to create GGt");
	createGGT_time.AddStart(omp_get_wtime());

	if (USE_DYNAMIC == 0) {
		if (USE_GGtINV == 1)
#ifdef DEVEL
			CreateGGt_inv_dist( cluster );
#else
			CreateGGt_inv( cluster );
#endif
		else
			CreateGGt    ( cluster );
	}

	createGGT_time.AddEnd(omp_get_wtime());
	createGGT_time.PrintStatMPI(0.0);
	preproc_timing.AddEvent(createGGT_time);



	// *** END - Make GGt *********************************************************
	// ****************************************************************************

	preproc_timing.totalTime.AddEnd(omp_get_wtime());


}

void IterSolver::Solve_singular ( Cluster & cluster, string & result_file )
{

#ifdef DEVEL
	if ( USE_PIPECG == 1 ) {
		Solve_PipeCG_singular_dom( cluster );
	} else {
		Solve_RegCG_singular_dom( cluster );
	}

	//MakeSolution_Primal_singular_parallel(cluster);
	//Save_to_Ensight_file( cluster, result_file );
#else
	if ( USE_PIPECG   == 1 )
		Solve_PipeCG_singular( cluster );
	else
		Solve_RegCG_singular ( cluster );

	if ( FIND_SOLUTION == 1 ) {
		MakeSolution_Primal_singular_parallel(cluster);
		Save_to_Ensight_file( cluster, result_file );
	}
#endif // DEVEL

}

void IterSolver::GetResiduum_Dual_singular_parallel    ( Cluster & cluster, SEQ_VECTOR <double> & dual_residuum_out ) {

	dual_residuum_out = dual_residuum_compressed_parallel;
	cluster.decompress_lambda_vector( dual_residuum_out );

}

void IterSolver::GetSolution_Dual_singular_parallel    ( Cluster & cluster, SEQ_VECTOR <double> & dual_solution_out, SEQ_VECTOR<double> & amplitudes_out ) {

	dual_solution_out = dual_soultion_compressed_parallel;
	cluster.decompress_lambda_vector( dual_solution_out );

	amplitudes_out	  = amplitudes;

}

void IterSolver::GetSolution_Primal_singular_parallel  ( Cluster & cluster, SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out ) {

	MakeSolution_Primal_singular_parallel(cluster);
	primal_solution_out = primal_solution_parallel;

}

// will be private
void IterSolver::MakeSolution_Primal_singular_parallel ( Cluster & cluster)  {

	primal_solution_parallel.clear();

	// R * mu
	SEQ_VECTOR<SEQ_VECTOR<double> > R_mu_prim_cluster;

	for (int d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR <double > tmp (cluster.domains[d].domain_prim_size);
		if (USE_HFETI == 1)
			cluster.domains[d].Kplus_R.MatVec(amplitudes, tmp, 'N', 0, 0);
		else
			cluster.domains[d].Kplus_R.MatVec(amplitudes, tmp, 'N', d * cluster.domains[d].Kplus_R.cols, 0);

		R_mu_prim_cluster.push_back(tmp);
	}


#ifdef DEVEL
	for (int d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1t_comp_dom.cols );
		SEQ_VECTOR < double > tmp      ( cluster.domains[d].domain_prim_size  );

		for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = dual_soultion_compressed_parallel[ cluster.domains[d].lambda_map_sub_local[i]];

		cluster.domains[d].B1t_comp_dom.MatVec (x_in_tmp, tmp, 'N');

		for (int i = 0; i < tmp.size(); i++)
			tmp[i] = cluster.domains[d].f[i] - tmp[i];

		primal_solution_parallel.push_back(tmp);

	}
#else
	for (int d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR <double > tmp (cluster.domains[d].domain_prim_size);
		cluster.domains[d].B1t_comp.MatVec(dual_soultion_compressed_parallel, tmp,'N');
		for (int i = 0; i < tmp.size(); i++)
			tmp[i] = cluster.domains[d].f[i] - tmp[i];
		primal_solution_parallel.push_back(tmp);
	}
#endif


	if ( cluster.USE_HFETI == 0) {
		for (int d = 0; d < cluster.domains.size(); d++)
			cluster.domains[d].multKplusLocal(primal_solution_parallel[d]);
	} else  {
		cluster.multKplusGlobal_l(primal_solution_parallel);
	}

 	for (int d = 0; d < cluster.domains.size(); d++) {
		for (int i = 0; i < primal_solution_parallel[d].size()	; i++) {
			primal_solution_parallel[d][i] = primal_solution_parallel[d][i] + R_mu_prim_cluster[d][i];
			//primal_solution_parallel[d][i] = cluster.domains[d].up0[i] + R_mu_prim_cluster[d][i];
		}
	}

	//// *** Saving results to binary file *******************************************
	//for (int d = 0; d < cluster.domains.size(); d++) {
	//	char number[8];
	//	sprintf(number, "%06d", cluster.domains[d].domain_global_index);
	//	string path = string(cluster.data_directory) + "/" + "res_u_" + number;
	//	SaveBinVectorDouble(primal_solution[d], path);
	//}
	//// *** END - Saving results to binary file *************************************


	//std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
	//SEQ_VECTOR < double > y_out_tmp (cluster.domains[0].B1_comp_dom.rows);

	//for (int d = 0; d < cluster.domains.size(); d++) {
	//	y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
	//	cluster.domains[d].B1_comp_dom.MatVec (primal_solution_parallel[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
	//	for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
	//		cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
	//}


	//vector < double > y_out (cluster.my_lamdas_indices.size(), 0.0);
	//All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	//double nnorm = parallel_norm_compressed(cluster,y_out);
	//
	//if (mpi_rank == 0)
	//	cout << "B*u norm = " << nnorm << endl;

}


// *** Output routines ***************************************************
void IterSolver::Save_to_Ensight_file ( Cluster & cluster, string & result_file) {

	int global_size = cluster.domains[0].number_of_nodes_in_global0[0];

	SEQ_VECTOR <double> global_ux (global_size, 0);
	SEQ_VECTOR <double> global_uy (global_size, 0);
	SEQ_VECTOR <double> global_uz (global_size, 0);

	for (int d = 0; d < cluster.domains.size(); d++) {
		for (int i = 0; i < primal_solution_parallel[d].size(); i = i + 3) {
			int node_index = i / 3;
			cluster.domains[d].ux.push_back( primal_solution_parallel[d][i  ] );	// ux
			cluster.domains[d].uy.push_back( primal_solution_parallel[d][i+1] );	// uy
			cluster.domains[d].uz.push_back( primal_solution_parallel[d][i+2] );    // uz

			global_ux[ cluster.domains[d].map_vector_local2global0[node_index] - 1 ] += primal_solution_parallel[d][i  ] / cluster.domains[d].nodeMulti[node_index];
			global_uy[ cluster.domains[d].map_vector_local2global0[node_index] - 1 ] += primal_solution_parallel[d][i+1] / cluster.domains[d].nodeMulti[node_index];
			global_uz[ cluster.domains[d].map_vector_local2global0[node_index] - 1 ] += primal_solution_parallel[d][i+2] / cluster.domains[d].nodeMulti[node_index];
		}
	}

	SEQ_VECTOR <double> final_global_ux (global_size, 0);
	SEQ_VECTOR <double> final_global_uy (global_size, 0);
	SEQ_VECTOR <double> final_global_uz (global_size, 0);

#ifdef LIB
	MPI_Reduce(&global_ux[0], &final_global_ux[0], global_size, MPI_DOUBLE, MPI_MAX, mpi_root, MPI_COMM_WORLD);
	MPI_Reduce(&global_uy[0], &final_global_uy[0], global_size, MPI_DOUBLE, MPI_MAX, mpi_root, MPI_COMM_WORLD);
	MPI_Reduce(&global_uz[0], &final_global_uz[0], global_size, MPI_DOUBLE, MPI_MAX, mpi_root, MPI_COMM_WORLD);
#else
	MPI_Reduce(&global_ux[0], &final_global_ux[0], global_size, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);
	MPI_Reduce(&global_uy[0], &final_global_uy[0], global_size, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);
	MPI_Reduce(&global_uz[0], &final_global_uz[0], global_size, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);
#endif
	if (mpi_rank == mpi_root) {

#ifdef WIN32
		_set_output_format(_TWO_DIGIT_EXPONENT);
#endif

		char   s0[] = "This is the description of the EnSight Gold geometry, MODEL:MATSOL_SVN";
		char   s1[] = "part";
		char   s2[] = "         1";
		char   s3[] = "coordinates";
		char   c = '\n';


		cout << "************************************************************************************************************************************************** " << endl;
		cout << " Saving results : " << endl;

		FILE *stream;
		stream = fopen( result_file.c_str(), "w" );

		fprintf( stream, "%s%c", s0, c );
		fprintf( stream, "%s%c", s1, c );
		fprintf( stream, "%s%c", s2, c );
		fprintf( stream, "%s%c", s3, c );

		for (int i = 0; i<final_global_ux.size(); i++)
			fprintf( stream, "%12.5e\n",final_global_ux[i]);

		for (int i = 0; i<final_global_ux.size(); i++)
			fprintf( stream, "%12.5e\n",final_global_uy[i]);

		for (int i = 0; i<final_global_ux.size(); i++)
			fprintf( stream, "%12.5e\n",final_global_uz[i]);

		fclose( stream );

		cout << "Output file is : " << result_file << endl;

	}

}

void IterSolver::Save_to_Ensight_file ( Cluster & cluster, string & result_file, SEQ_VECTOR < SEQ_VECTOR < double> > & in_primal_solution_parallel ) {

	int global_size = cluster.domains[0].number_of_nodes_in_global0[0];
	SEQ_VECTOR <double> global_ux (global_size, 0);
	SEQ_VECTOR <double> global_uy (global_size, 0);
	SEQ_VECTOR <double> global_uz (global_size, 0);

	for (int d = 0; d < cluster.domains.size(); d++) {
		for (int i = 0; i < in_primal_solution_parallel[d].size(); i = i + 3) {
			int node_index = i / 3;

			cluster.domains[d].ux.push_back( in_primal_solution_parallel[d][i  ] );	// ux
			cluster.domains[d].uy.push_back( in_primal_solution_parallel[d][i+1] );	// uy
			cluster.domains[d].uz.push_back( in_primal_solution_parallel[d][i+2] );	// uz

			global_ux[ cluster.domains[d].map_vector_local2global0[node_index] - 1 ] += in_primal_solution_parallel[d][i  ] / cluster.domains[d].nodeMulti[node_index];
			global_uy[ cluster.domains[d].map_vector_local2global0[node_index] - 1 ] += in_primal_solution_parallel[d][i+1] / cluster.domains[d].nodeMulti[node_index];
			global_uz[ cluster.domains[d].map_vector_local2global0[node_index] - 1 ] += in_primal_solution_parallel[d][i+2] / cluster.domains[d].nodeMulti[node_index];


		}
	}

	SEQ_VECTOR <double> final_global_ux (global_size, 0);
	SEQ_VECTOR <double> final_global_uy (global_size, 0);
	SEQ_VECTOR <double> final_global_uz (global_size, 0);

	MPI_Reduce(&global_ux[0], &final_global_ux[0], global_size, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);
	MPI_Reduce(&global_uy[0], &final_global_uy[0], global_size, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);
	MPI_Reduce(&global_uz[0], &final_global_uz[0], global_size, MPI_DOUBLE, MPI_SUM, mpi_root, MPI_COMM_WORLD);

	//save to file
	if (mpi_rank == mpi_root) {

#ifdef WIN32
		_set_output_format(_TWO_DIGIT_EXPONENT);
#endif

		char   s0[] = "This is the description of the EnSight Gold geometry, MODEL:MATSOL_SVN";
		char   s1[] = "part";
		char   s2[] = "         1";
		char   s3[] = "coordinates";
		char   c = '\n';

		FILE *stream;
		stream = fopen( result_file.c_str(), "w" );

		fprintf( stream, "%s%c", s0, c );
		fprintf( stream, "%s%c", s1, c );
		fprintf( stream, "%s%c", s2, c );
		fprintf( stream, "%s%c", s3, c );

		for (int i = 0; i<final_global_ux.size(); i++)
			fprintf( stream, "%12.5e\n",final_global_ux[i]);

		for (int i = 0; i<final_global_ux.size(); i++)
			fprintf( stream, "%12.5e\n",final_global_uy[i]);

		for (int i = 0; i<final_global_ux.size(); i++)
			fprintf( stream, "%12.5e\n",final_global_uz[i]);

		fclose( stream );
	}

}


// *** Singular CG Solvers ***********************************************
void IterSolver::Solve_RegCG_singular ( Cluster & cluster ) //, vector <double> & x_l)   // dual_soultion_in = x_l
{

	int dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);

	SEQ_VECTOR <double> Ax_l (dl_size, 0);
	SEQ_VECTOR <double> Ap_l (dl_size, 0);
	SEQ_VECTOR <double> r_l  (dl_size, 0);

	SEQ_VECTOR <double> w_l  (dl_size, 0);
	SEQ_VECTOR <double> wp_l (dl_size, 0);

	SEQ_VECTOR <double> y_l  (dl_size, 0);
	SEQ_VECTOR <double> yp_l (dl_size, 0);
	SEQ_VECTOR <double> z_l  (dl_size, 0);
	SEQ_VECTOR <double> p_l  (dl_size, 0);

	SEQ_VECTOR <double> u_l  (dl_size, 0);

	SEQ_VECTOR <double> b_l  (dl_size, 0);

	double beta_l  = 0;
	double alpha_l = 0;
	double norm_l;
	double tol;

	// *** CG start ***************************************************************

	// t1 = Uc\(Lc\d);
	// x = Ct * t1;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	// *** r = b - Ax *************************************************************
	cilk_for (int i = 0; i < r_l.size(); i++)
		r_l[i] = b_l[i] - Ax_l[i];

	if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, r_l, u_l , 0);
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, r_l, u_l , 0);
		}

	// *** Calculate the stop condition *******************************************
	tol = epsilon * parallel_norm_compressed(cluster, u_l);

	// *** Start the CG iteration loop ********************************************
	for (int iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.AddStart(omp_get_wtime());

		cilk_for (int i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		if (USE_PREC == 1) {

			proj1_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj1_time.AddEnd(omp_get_wtime());


			prec_time.AddStart(omp_get_wtime());
			apply_prec_compB(timeEvalPrec, cluster, w_l, z_l);
			prec_time.AddEnd(omp_get_wtime());


			proj2_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, z_l, y_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, z_l, y_l, 0 );
			}
			proj2_time.AddEnd(omp_get_wtime());

		} else {

			proj_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj_time.AddEnd(omp_get_wtime());

			cilk_for (int i = 0; i < w_l.size(); i++)
				y_l[i] = w_l[i];

		}


		//------------------------------------------
		if (iter == 0) {									// if outputs.n_it==1;

			cilk_for (int i = 0; i < y_l.size(); i++)
				p_l[i] = y_l[i];							// p = y;

		} else {

			ddot_beta.AddStart(omp_get_wtime());
			beta_l =          parallel_ddot_compressed(cluster, y_l, w_l);
			beta_l = beta_l / parallel_ddot_compressed(cluster, yp_l, wp_l);
			ddot_beta.AddEnd(omp_get_wtime());

			cilk_for (int i = 0; i < p_l.size(); i++)
				p_l[i] = y_l[i] + beta_l * p_l[i];			// p = y + beta * p;

		}



		//------------------------------------------
		appA_time.AddStart(omp_get_wtime());
		apply_A_l_compB(timeEvalAppa, cluster, p_l, Ap_l);
		appA_time.AddEnd(omp_get_wtime());

		//------------------------------------------
		ddot_alpha.AddStart(omp_get_wtime());
		alpha_l = parallel_ddot_compressed(cluster, y_l, w_l);
		alpha_l = alpha_l / parallel_ddot_compressed(cluster, p_l, Ap_l);
		ddot_alpha.AddEnd(omp_get_wtime());

		//------------------------------------------
		cilk_for (int i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * Ap_l[i];
		}

		timing.totalTime.AddEnd(omp_get_wtime());
		timing.totalTime.PrintLastStatMPI(0.0);

		norm_time.AddStart(omp_get_wtime());
		norm_l = parallel_norm_compressed(cluster, w_l);
		norm_time.AddEnd(omp_get_wtime());

		if (mpi_rank == mpi_root) {
			printf (       "Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);

			//if (log_active == 1)
				//fprintf(stream,"Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
		}


		// *** Stop condition ******************************************************************
		if (norm_l < tol)
			break;

	} // end of CG iterations


	// *** save solution - in dual and amplitudes *********************************************
	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	//if (mpi_rank == mpi_root)
		//if (log_active == 1)
		//	fclose( stream );

	// *** Print out the timing for the iteration loop ***************************************

	if (USE_PREC == 1) {
		timing.AddEvent(proj1_time);
		timing.AddEvent(prec_time );
		timing.AddEvent(proj2_time);
	} else {
		timing.AddEvent(proj_time);
	}

	timing.AddEvent(appA_time );
	timing.AddEvent(ddot_beta);
	timing.AddEvent(ddot_alpha);

	preproc_timing.PrintStatsMPI();
	timing.PrintStatsMPI();
	timeEvalAppa.PrintStatsMPI();
	timeEvalProj.PrintStatsMPI();

	if (USE_PREC == 1)
		timeEvalPrec.PrintStatsMPI();

	//if(cluster.domains.size() > 1 )
	if ( cluster.USE_HFETI == 1 )
		cluster.ShowTiming();

	// *** END - Print out the timing for the iteration loop ***********************************

}

void IterSolver::Solve_RegCG_singular_dom ( Cluster & cluster ) //, vector <double> & x_l)   // dual_soultion_in = x_l
{

	int dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);

	SEQ_VECTOR <double> Ax_l (dl_size, 0);
	SEQ_VECTOR <double> Ap_l (dl_size, 0);
	SEQ_VECTOR <double> r_l  (dl_size, 0);

	SEQ_VECTOR <double> w_l  (dl_size, 0);
	SEQ_VECTOR <double> wp_l (dl_size, 0);

	SEQ_VECTOR <double> y_l  (dl_size, 0);
	SEQ_VECTOR <double> yp_l (dl_size, 0);
	SEQ_VECTOR <double> z_l  (dl_size, 0);
	SEQ_VECTOR <double> p_l  (dl_size, 0);

	SEQ_VECTOR <double> u_l  (dl_size, 0);

	SEQ_VECTOR <double> b_l  (dl_size, 0);

	double beta_l  = 0;
	double alpha_l = 0;
	double norm_l;
	double tol;

	// *** CG start ***************************************************************

	// t1 = Uc\(Lc\d);
	// x = Ct * t1;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	// *** up0 pro ukoncovani v primaru
	////cilk_for (int d = 0; d < cluster.domains.size(); d++) {
	////	cluster.domains[d].BtLambda_i = cluster.x_prim_cluster2[d];
	////	//cluster.domains[d].BtLambda_i.resize(cluster.domains[d].up0.size(), 0);
	////	cluster.domains[d].norm_f = 0.0;
	////	for (int i = 0; i < cluster.domains[d].up0.size(); i++ ) {
	////		cluster.domains[d].up0[i]     = cluster.domains[d].up0[i] - cluster.x_prim_cluster1[d][i];  // (K+ * f) - (K+ * Bt * lambda)
	////
	////		cluster.domains[d].norm_f += cluster.domains[d].f[i] * cluster.domains[d].f[i];
	////	}
	////}

	// Get norm of f (right hand side)
	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (int d = 0; d < cluster.domains.size(); d++)
		norm_prim_fl += cluster.domains[d].norm_f;

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Reduce   (&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	norm_prim_fg = sqrt(norm_prim_fg);



	// *** r = b - Ax *************************************************************
	cilk_for (int i = 0; i < r_l.size(); i++)
		r_l[i] = b_l[i] - Ax_l[i];

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, r_l, u_l , 0);
	} else {
		Projector_l_compG    ( timeEvalProj, cluster, r_l, u_l , 0);
	}

	// *** Calculate the stop condition *******************************************
	tol = epsilon * parallel_norm_compressed(cluster, u_l);

	// *** Start the CG iteration loop ********************************************
	for (int iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.AddStart(omp_get_wtime());

		cilk_for (int i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		if (USE_PREC == 1) {

			proj1_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj1_time.AddEnd(omp_get_wtime());

			// Scale
			prec_time.AddStart(omp_get_wtime());
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, z_l);
			prec_time.AddEnd(omp_get_wtime());
			// Re-Scale

			proj2_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, z_l, y_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, z_l, y_l, 0 );
			}
			proj2_time.AddEnd(omp_get_wtime());

		} else {

			proj_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj_time.AddEnd(omp_get_wtime());

			cilk_for (int i = 0; i < w_l.size(); i++)
				y_l[i] = w_l[i];

		}


		//------------------------------------------
		if (iter == 0) {									// if outputs.n_it==1;

			cilk_for (int i = 0; i < y_l.size(); i++)
				p_l[i] = y_l[i];							// p = y;

		} else {

			ddot_beta.AddStart(omp_get_wtime());
			beta_l =          parallel_ddot_compressed(cluster, y_l, w_l);
			beta_l = beta_l / parallel_ddot_compressed(cluster, yp_l, wp_l);
			ddot_beta.AddEnd(omp_get_wtime());

			cilk_for (int i = 0; i < p_l.size(); i++)
				p_l[i] = y_l[i] + beta_l * p_l[i];			// p = y + beta * p;

		}



		//------------------------------------------
		 appA_time.AddStart(omp_get_wtime());
		apply_A_l_comp_dom_B(timeEvalAppa, cluster, p_l, Ap_l); // apply_A_l_compB(timeEvalAppa, cluster, p_l, Ap_l);
		 appA_time.AddEnd(omp_get_wtime());

		//------------------------------------------
		 ddot_alpha.AddStart(omp_get_wtime());
		alpha_l =           parallel_ddot_compressed(cluster, y_l, w_l);
		alpha_l = alpha_l / parallel_ddot_compressed(cluster, p_l, Ap_l);
		 ddot_alpha.AddEnd(omp_get_wtime());

		//-----------------------------------------
		// *** up0 pro ukoncovani v primaru

		//// //cilk_
		////for (int d = 0; d < cluster.domains.size(); d++) {
		////	for (int i = 0; i < cluster.domains[d].up0.size(); i++ ) {
		////		cluster.domains[d].up0[i]        -= alpha_l * cluster.x_prim_cluster1[d][i];
		////		cluster.domains[d].BtLambda_i[i] += alpha_l * cluster.x_prim_cluster2[d][i];
		////	}
		////	cluster.domains[d].norm_vec.resize(cluster.domains[d].up0.size());
		////	cluster.domains[d].K.MatVec(cluster.domains[d].up0, cluster.domains[d].norm_vec, 'N');
		////	cluster.domains[d].norm_c = 0.0;
		////	for (int i = 0; i < cluster.domains[d].up0.size(); i++ ) {
		////		cluster.domains[d].norm_vec[i] = cluster.domains[d].norm_vec[i]
		////			                           + cluster.domains[d].BtLambda_i[i]
		////									   - cluster.domains[d].f[i];
		////
		////		cluster.domains[d].norm_c += cluster.domains[d].norm_vec[i] * cluster.domains[d].norm_vec[i];
		////	}
		////}

		//double norm_prim_l = 0.0;
		//double norm_prim_g = 0.0;
		//for (int d = 0; d < cluster.domains.size(); d++)
		//	norm_prim_l += cluster.domains[d].norm_c;

		//MPI_Allreduce(&norm_prim_l, &norm_prim_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		////MPI_Reduce(&norm_prim_l, &norm_prim_g, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		//norm_prim_g = sqrt(norm_prim_g);





		//------------------------------------------
		cilk_for (int i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * Ap_l[i];
		}

		timing.totalTime.AddEnd(omp_get_wtime());
		timing.totalTime.PrintLastStatMPI(0.0);

		 norm_time.AddStart(omp_get_wtime());
		norm_l = parallel_norm_compressed(cluster, w_l);
		 norm_time.AddEnd(omp_get_wtime());


		if (mpi_rank == mpi_root) {
			//printf (       "Iter MPI %d - norm dual %1.20f - tol %1.20f - norm prim %1.20f norm f %1.20f : %1.20f \n", iter+1, norm_l, tol, norm_prim_g, norm_prim_fg, norm_prim_g / norm_prim_fg);
			printf (       "Iter MPI %d - norm dual %1.20f - tol %1.20f \n", iter+1, norm_l, tol );

			//if (log_active == 1)
			//fprintf(stream,"Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
		}


		// *** Stop condition ******************************************************************
		if (norm_l < tol)
			break;

	} // end of CG iterations


	// *** save solution - in dual and amplitudes *********************************************
	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************

	//MakeSolution_Primal_singular_parallel(cluster);


	//if (mpi_rank == mpi_root)
	//if (log_active == 1)
	//	fclose( stream );

	// *** Print out the timing for the iteration loop ***************************************

	if (USE_PREC == 1) {
		timing.AddEvent(proj1_time);
		timing.AddEvent(prec_time );
		timing.AddEvent(proj2_time);
	} else {
		timing.AddEvent(proj_time);
	}

	timing.AddEvent(appA_time );
	timing.AddEvent(ddot_beta);
	timing.AddEvent(ddot_alpha);

	preproc_timing.PrintStatsMPI();
	timing.PrintStatsMPI();
	timeEvalAppa.PrintStatsMPI();
	timeEvalProj.PrintStatsMPI();

	if (USE_PREC == 1)
		timeEvalPrec.PrintStatsMPI();

	//if(cluster.domains.size() > 1 )
	if ( cluster.USE_HFETI == 1 )
		cluster.ShowTiming();

	// *** END - Print out the timing for the iteration loop ***********************************

}


void IterSolver::Solve_PipeCG_singular ( Cluster & cluster ) //, vector <double> & x_l) // dual_soultion_in = x_l
{
	int dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);

	SEQ_VECTOR <double> b_l  (dl_size, 0);

	SEQ_VECTOR <double> Ax_l (dl_size, 0);
	SEQ_VECTOR <double> Ap_l (dl_size, 0);
	SEQ_VECTOR <double> r_l  (dl_size, 0);

	SEQ_VECTOR <double> w_l  (dl_size, 0);
	SEQ_VECTOR <double> wp_l (dl_size, 0);

	SEQ_VECTOR <double> y_l  (dl_size, 0);
	SEQ_VECTOR <double> yp_l (dl_size, 0);

	SEQ_VECTOR <double> z_l  (dl_size, 0);
	SEQ_VECTOR <double> p_l  (dl_size, 0);

	SEQ_VECTOR <double> s_l  (dl_size, 0);
	SEQ_VECTOR <double> q_l  (dl_size, 0);
	SEQ_VECTOR <double> m_l  (dl_size, 0);
	SEQ_VECTOR <double> n_l  (dl_size, 0);

	SEQ_VECTOR <double> u_l  (dl_size, 0);
	SEQ_VECTOR <double> tmp_l(dl_size, 0);


	double gama_l  = 0;
	double gama_lp = 0;

	double delta_l = 0;
	double beta_l  = 0;

	double alpha_l = 0;
	double alpha_lp;

	double norm_l;
	double tol;

	// *** CG start ***************************************************************
	// t1 = Uc\(Lc\d);
	// x = Ct * t1;

	if (USE_GGtINV == 1)
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	else
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );


	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);


	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);


	// *** r = b - Ax; ************************************************************
	if (USE_PREC == 1) {

		cilk_for (int i = 0; i < r_l.size(); i++)
			tmp_l[i] = b_l[i] - Ax_l[i];

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, r_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, r_l, 0 );

	} else {

		cilk_for (int i = 0; i < r_l.size(); i++)
			r_l[i] = b_l[i] - Ax_l[i];
	}


	if (USE_PREC == 1) {

		apply_prec_compB(timeEvalPrec, cluster, r_l, tmp_l);

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, u_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, u_l, 0 );

	} else {

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, r_l, u_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, r_l, u_l, 0 );
	}


	if (USE_PREC == 1) {
		apply_A_l_compB(timeEvalAppa, cluster, u_l, tmp_l);

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, w_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, w_l, 0 );

	} else {

		apply_A_l_compB(timeEvalAppa, cluster, u_l, w_l);

	}


	// *** Calculate the stop condition *******************************************

	if (USE_GGtINV == 1)
		Projector_l_inv_compG( timeEvalProj, cluster, b_l, tmp_l, 0 );
	else
		Projector_l_compG    ( timeEvalProj, cluster, b_l, tmp_l, 0 );

	tol = epsilon * parallel_norm_compressed(cluster, tmp_l);
	//tol = epsilon * parallel_norm_compressed(cluster, u_l);



	// *** Start the CG iteration loop ********************************************
	for (int iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.AddStart(omp_get_wtime());

		alpha_lp = alpha_l;
		gama_lp  = gama_l;

		//------------------------------------------
		ddot_time.AddStart(omp_get_wtime());
		MPI_Request mpi_req;
		MPI_Status mpi_stat;

		SEQ_VECTOR <double> reduction_tmp (2,0);
		SEQ_VECTOR <double> send_buf (2,0);
		parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, &mpi_req, reduction_tmp, send_buf);

		ddot_time.AddEnd(omp_get_wtime());

		if (USE_PREC == 1) {

			prec_time.AddStart(omp_get_wtime());
			apply_prec_compB(timeEvalPrec, cluster, w_l, tmp_l);
			prec_time.AddEnd(omp_get_wtime());

			proj_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, m_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, m_l, 0 );
			proj_time.AddEnd(omp_get_wtime());

			appA_time.AddStart(omp_get_wtime());
			apply_A_l_compB(timeEvalAppa, cluster, m_l, tmp_l);
			appA_time.AddEnd(omp_get_wtime());

			proj_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, n_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, n_l, 0 );
			proj_time.AddEnd(omp_get_wtime());

		} else {

			//------------------------------------------
			proj_time.AddStart(omp_get_wtime());

			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, w_l, m_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, w_l, m_l, 0 );

			proj_time.AddEnd(omp_get_wtime());

			//------------------------------------------
			appA_time.AddStart(omp_get_wtime());
			apply_A_l_compB(timeEvalAppa, cluster, m_l, n_l);
			appA_time.AddEnd(omp_get_wtime());
			//------------------------------------------
		}


#ifndef WIN32
#ifdef USE_MPI_3
		MPI_Wait(&mpi_req, &mpi_stat);
#endif
#endif
		gama_l  = reduction_tmp[0];
		delta_l = reduction_tmp[1];

		//------------------------------------------
		vec_time.AddStart(omp_get_wtime());
		if (iter == 0) {
			beta_l  = 0;
			alpha_l = gama_l / delta_l;
		} else {
			beta_l  = gama_l / gama_lp;
			alpha_l = gama_l / (delta_l - beta_l * gama_l / alpha_lp);
		}

		cilk_for (int i = 0; i < r_l.size(); i++) {
			z_l[i] = n_l[i] + beta_l  * z_l[i];
			q_l[i] = m_l[i] + beta_l  * q_l[i];
			s_l[i] = w_l[i] + beta_l  * s_l[i];
			p_l[i] = u_l[i] + beta_l  * p_l[i];
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * s_l[i];
			u_l[i] = u_l[i] - alpha_l * q_l[i];
			w_l[i] = w_l[i] - alpha_l * z_l[i];
		}
		vec_time.AddEnd(omp_get_wtime());


		timing.totalTime.AddEnd(omp_get_wtime());
		timing.totalTime.PrintLastStatMPI(0.0);

		norm_time.AddStart(omp_get_wtime());

		// POZOR - tady se to ukoncuje jinak = musime probrat
		if (USE_PREC == 1)
			norm_l = parallel_norm_compressed(cluster, r_l);
		else
			norm_l = parallel_norm_compressed(cluster, u_l);

		norm_time.AddEnd(omp_get_wtime());

		if (mpi_rank == mpi_root) {
			printf (       "Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);

			//if (log_active == 1)
			//fprintf(stream,"Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
		}

		if (norm_l < tol)
			break;
	}

	// *** save solution - in dual and amplitudes *********************************************

	apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	cilk_for(int i = 0; i < r_l.size(); i++)
		r_l[i] = b_l[i] - Ax_l[i];

	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	//if (mpi_rank == mpi_root)
	//if (log_active == 1)
	//	fclose( stream );

	// *** Print out the timing for the iteration loop ***************************************

	timing.AddEvent(ddot_time);
	timing.AddEvent(proj_time);
	timing.AddEvent(appA_time);
	timing.AddEvent(vec_time );

	preproc_timing.PrintStatsMPI();
	timing.PrintStatsMPI();
	timeEvalAppa.PrintStatsMPI();
	timeEvalProj.PrintStatsMPI();

	if (USE_PREC == 1)
		timeEvalPrec.PrintStatsMPI();

	//if(cluster.domains.size() > 1 )
	if ( cluster.USE_HFETI == 1 )
		cluster.ShowTiming();


}

void IterSolver::Solve_PipeCG_singular_dom ( Cluster & cluster ) //, vector <double> & x_l) // dual_soultion_in = x_l
{
	int dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);

	SEQ_VECTOR <double> b_l  (dl_size, 0);

	SEQ_VECTOR <double> Ax_l (dl_size, 0);
	SEQ_VECTOR <double> Ap_l (dl_size, 0);
	SEQ_VECTOR <double> r_l  (dl_size, 0);

	SEQ_VECTOR <double> w_l  (dl_size, 0);
	SEQ_VECTOR <double> wp_l (dl_size, 0);

	SEQ_VECTOR <double> y_l  (dl_size, 0);
	SEQ_VECTOR <double> yp_l (dl_size, 0);

	SEQ_VECTOR <double> z_l  (dl_size, 0);
	SEQ_VECTOR <double> p_l  (dl_size, 0);

	SEQ_VECTOR <double> s_l  (dl_size, 0);
	SEQ_VECTOR <double> q_l  (dl_size, 0);
	SEQ_VECTOR <double> m_l  (dl_size, 0);
	SEQ_VECTOR <double> n_l  (dl_size, 0);

	SEQ_VECTOR <double> u_l  (dl_size, 0);
	SEQ_VECTOR <double> tmp_l(dl_size, 0);


	double gama_l  = 0;
	double gama_lp = 0;

	double delta_l = 0;
	double beta_l  = 0;

	double alpha_l = 0;
	double alpha_lp;

	double norm_l;
	double tol;

	// *** CG start ***************************************************************
	// t1 = Uc\(Lc\d);
	// x = Ct * t1;

	if (USE_GGtINV == 1)
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	else
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );


	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);


	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l); //apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);


	// *** r = b - Ax; ************************************************************
	if (USE_PREC == 1) {

		cilk_for (int i = 0; i < r_l.size(); i++)
			tmp_l[i] = b_l[i] - Ax_l[i];

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, r_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, r_l, 0 );

	} else {

		cilk_for (int i = 0; i < r_l.size(); i++)
			r_l[i] = b_l[i] - Ax_l[i];
	}


	if (USE_PREC == 1) {

		apply_prec_comp_dom_B(timeEvalPrec, cluster, r_l, tmp_l);

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, u_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, u_l, 0 );

	} else {

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, r_l, u_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, r_l, u_l, 0 );
	}


	if (USE_PREC == 1) {
		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, tmp_l); //apply_A_l_compB(timeEvalAppa, cluster, u_l, tmp_l);
		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, w_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, w_l, 0 );

	} else {

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, w_l); 	//apply_A_l_compB(timeEvalAppa, cluster, u_l, w_l);

	}


	// *** Calculate the stop condition *******************************************

	if (USE_GGtINV == 1)
		Projector_l_inv_compG( timeEvalProj, cluster, b_l, tmp_l, 0 );
	else
		Projector_l_compG    ( timeEvalProj, cluster, b_l, tmp_l, 0 );

	tol = epsilon * parallel_norm_compressed(cluster, tmp_l); //tol = epsilon * parallel_norm_compressed(cluster, u_l);



	// *** Start the CG iteration loop ********************************************
	for (int iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.AddStart(omp_get_wtime());

		alpha_lp = alpha_l;
		gama_lp  = gama_l;

		//------------------------------------------
		ddot_time.AddStart(omp_get_wtime());
		MPI_Request mpi_req;
		MPI_Status mpi_stat;

		SEQ_VECTOR <double> reduction_tmp (2,0);
		SEQ_VECTOR <double> send_buf (2,0);
		parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, &mpi_req, reduction_tmp, send_buf);

		ddot_time.AddEnd(omp_get_wtime());

		if (USE_PREC == 1) {

			prec_time.AddStart(omp_get_wtime());
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, tmp_l);
			prec_time.AddEnd(omp_get_wtime());

			proj_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, m_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, m_l, 0 );
			proj_time.AddEnd(omp_get_wtime());

			appA_time.AddStart(omp_get_wtime());
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, tmp_l); //apply_A_l_compB(timeEvalAppa, cluster, m_l, tmp_l);
			appA_time.AddEnd(omp_get_wtime());

			proj_time.AddStart(omp_get_wtime());
			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, n_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, n_l, 0 );
			proj_time.AddEnd(omp_get_wtime());

		} else {

			//------------------------------------------
			proj_time.AddStart(omp_get_wtime());

			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, w_l, m_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, w_l, m_l, 0 );

			proj_time.AddEnd(omp_get_wtime());

			//------------------------------------------
			appA_time.AddStart(omp_get_wtime());
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, n_l); //apply_A_l_compB(timeEvalAppa, cluster, m_l, n_l);
			appA_time.AddEnd(omp_get_wtime());
			//------------------------------------------
		}

#ifndef WIN32
#ifdef USE_MPI_3
		MPI_Wait(&mpi_req, &mpi_stat);
#endif
#endif
		gama_l  = reduction_tmp[0];
		delta_l = reduction_tmp[1];

		//------------------------------------------
		vec_time.AddStart(omp_get_wtime());
		if (iter == 0) {
			beta_l  = 0;
			alpha_l = gama_l / delta_l;
		} else {
			beta_l  = gama_l / gama_lp;
			alpha_l = gama_l / (delta_l - beta_l * gama_l / alpha_lp);
		}

		cilk_for (int i = 0; i < r_l.size(); i++) {
			z_l[i] = n_l[i] + beta_l  * z_l[i];
			q_l[i] = m_l[i] + beta_l  * q_l[i];
			s_l[i] = w_l[i] + beta_l  * s_l[i];
			p_l[i] = u_l[i] + beta_l  * p_l[i];
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * s_l[i];
			u_l[i] = u_l[i] - alpha_l * q_l[i];
			w_l[i] = w_l[i] - alpha_l * z_l[i];
		}
		vec_time.AddEnd(omp_get_wtime());


		timing.totalTime.AddEnd(omp_get_wtime());
		timing.totalTime.PrintLastStatMPI(0.0);

		norm_time.AddStart(omp_get_wtime());

		// POZOR - tady se to ukoncuje jinak = musime probrat
		if (USE_PREC == 1)
			norm_l = parallel_norm_compressed(cluster, r_l);
		else
			norm_l = parallel_norm_compressed(cluster, u_l);

		norm_time.AddEnd(omp_get_wtime());

		if (mpi_rank == mpi_root) {
			printf (       "Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);

			//if (log_active == 1)
			//fprintf(stream,"Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
		}

		if (norm_l < tol)
			break;
	}

	// *** save solution - in dual and amplitudes *********************************************

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l); //apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	cilk_for(int i = 0; i < r_l.size(); i++)
		r_l[i] = b_l[i] - Ax_l[i];

	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	//MakeSolution_Primal_singular_parallel(cluster);

	//if (mpi_rank == mpi_root)
	//if (log_active == 1)
	//	fclose( stream );

	// *** Print out the timing for the iteration loop ***************************************

	timing.AddEvent(ddot_time);
	timing.AddEvent(proj_time);
	timing.AddEvent(appA_time);
	timing.AddEvent(vec_time );

	preproc_timing.PrintStatsMPI();
	timing.PrintStatsMPI();
	timeEvalAppa.PrintStatsMPI();
	timeEvalProj.PrintStatsMPI();

	if (USE_PREC == 1)
		timeEvalPrec.PrintStatsMPI();

	//if(cluster.domains.size() > 1 )
	if ( cluster.USE_HFETI == 1 )
		cluster.ShowTiming();


}


// *** Non-singular CG Solvers *******************************************
void IterSolver::Solve_RegCG_nonsingular  ( Cluster & cluster,
										    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
										    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel) {

	int dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l  (dl_size, 0);
	SEQ_VECTOR <double> b_l  (dl_size, 0);

	SEQ_VECTOR <double> Ax_l (dl_size, 0);
	SEQ_VECTOR <double> Ap_l (dl_size, 0);
	SEQ_VECTOR <double> r_l  (dl_size, 0);

	SEQ_VECTOR <double> w_l  (dl_size, 0);
	SEQ_VECTOR <double> wp_l (dl_size, 0);

	SEQ_VECTOR <double> y_l  (dl_size, 0);
	SEQ_VECTOR <double> yp_l (dl_size, 0);
	SEQ_VECTOR <double> z_l  (dl_size, 0);
	SEQ_VECTOR <double> p_l  (dl_size, 0);

	SEQ_VECTOR <double> s_l  (dl_size, 0);
	SEQ_VECTOR <double> q_l  (dl_size, 0);
	SEQ_VECTOR <double> m_l  (dl_size, 0);
	SEQ_VECTOR <double> n_l  (dl_size, 0);

	SEQ_VECTOR <double> u_l  (dl_size, 0);

	double gama_l  = 0;
	double gama_lp = 0;

	double delta_l = 0;
	double beta_l  = 0;

	double alpha_l = 0;
	double alpha_lp;

	double norm_l;
	double tol;


    // *** CG iteration loop ******************************************************************
	// *** - input to CG - vec_b - right hand side in primal


//	cluster.domains[0].B1_comp.MatVec(cluster.x_prim_cluster1[0], cluster.compressed_tmp, 'N', 0, 0, 0.0);
//	for (int d = 1; d < cluster.domains.size(); d++) {
//		cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[0]);
//		cluster.domains[d].B1_comp.MatVec(cluster.x_prim_cluster1[0], cluster.compressed_tmp, 'N', 0, 0, 1.0);
//	}
//	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);

	//// *** convert right hand side to dual
	std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
	SEQ_VECTOR < double > y_out_tmp;
	for (int d = 0; d < cluster.domains.size(); d++) {
		// *** convert right hand side to dual
		cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[d]);

		y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
		cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

		for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
	}
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);



	int iter = 0;
	timing.totalTime.Reset();

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);

	cilk_for (int i = 0; i < r_l.size(); i++) {	// r = b - Ax;
		r_l[i] = b_l[i] - Ax_l[i];
		wp_l[i] = 0.0;
		yp_l[i] = 0.0;
	}

	tol = epsilon * parallel_norm_compressed(cluster, b_l);

	for (iter = 0; iter < 1000; iter++) {
		timing.totalTime.AddStart(omp_get_wtime());

		cilk_for (int i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		if (USE_PREC == 1) {
			cilk_for (int i = 0; i < w_l.size(); i++)
				w_l[i] = r_l[i];
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, y_l);
		} else {
			cilk_for (int i = 0; i < w_l.size(); i++)
				w_l[i] = r_l[i];
			cilk_for (int i = 0; i < w_l.size(); i++)
				y_l[i] = w_l[i];
		}


		if (iter == 0) {									// if outputs.n_it==1;
			cilk_for (int i = 0; i < y_l.size(); i++)
				p_l[i] = y_l[i];							// p = y;
		} else {
			ddot_beta.AddStart(omp_get_wtime());
			beta_l =          parallel_ddot_compressed(cluster, y_l, w_l);
			beta_l = beta_l / parallel_ddot_compressed(cluster, yp_l, wp_l);
			ddot_beta.AddEnd(omp_get_wtime());

			cilk_for (int i = 0; i < p_l.size(); i++)
				p_l[i] = y_l[i] + beta_l * p_l[i];			// p = y + beta * p;
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, p_l, Ap_l);



		alpha_l =           parallel_ddot_compressed(cluster, y_l, w_l);
		alpha_l = alpha_l / parallel_ddot_compressed(cluster, p_l, Ap_l);

		cilk_for (int i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * Ap_l[i];
		}

		timing.totalTime.AddEnd(omp_get_wtime());

		norm_l = parallel_norm_compressed(cluster, r_l);

		if (mpi_rank == mpi_root) {
			printf ( "Iter MPI %d - norm %1.30f - tol %1.30f \r" , iter+1, norm_l, tol);
			//if (log_active == 1)
			//	fprintf(stream,"Iter MPI %d - norm %1.30f - tol %1.30f \n", iter+1, norm_l, tol);
		}


		if (norm_l < tol)
			break;

	} // end iter loop


	//// reconstruction of u
	//cluster.domains[0].B1t_comp.MatVec(x_l, vec_b_tmp, 'N');

	//cilk_for(int i = 0; i < in_right_hand_side_primal.size(); i++)
	//	out_primal_solution_parallel[i] = in_right_hand_side_primal[i] - vec_b_tmp[i];

	//cluster.domains[0].multKplusLocal(out_primal_solution_parallel);

	//// *** output of CG - vec_u_n in primal


	// reconstruction of u
	//cilk_
	for (int d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1t_comp_dom.cols, 0.0 );

		for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_l[ cluster.domains[d].lambda_map_sub_local[i]];

		cluster.domains[d].B1t_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');

		for(int i = 0; i < in_right_hand_side_primal[d].size(); i++)
			out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];

		cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);

	}

	// -- old --
//	for (int d = 0; d < cluster.domains.size(); d++) {
//		cluster.domains[d].B1t_comp.MatVec(x_l, cluster.x_prim_cluster1[d], 'N');
//
//		cilk_for(int i = 0; i < in_right_hand_side_primal[d].size(); i++)
//			out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];
//
//		cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);
//	}
	// *** output of CG - vec_u_n in primal


}

void IterSolver::Solve_PipeCG_nonsingular ( Cluster & cluster,
											SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
											SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel) {

		int dl_size = cluster.my_lamdas_indices.size();

		SEQ_VECTOR <double> x_l  (dl_size, 0);
		SEQ_VECTOR <double> b_l  (dl_size, 0);

		SEQ_VECTOR <double> Ax_l (dl_size, 0);
		SEQ_VECTOR <double> Ap_l (dl_size, 0);
		SEQ_VECTOR <double> r_l  (dl_size, 0);

		SEQ_VECTOR <double> w_l  (dl_size, 0);
		SEQ_VECTOR <double> wp_l (dl_size, 0);

		SEQ_VECTOR <double> y_l  (dl_size, 0);
		SEQ_VECTOR <double> yp_l (dl_size, 0);
		SEQ_VECTOR <double> z_l  (dl_size, 0);
		SEQ_VECTOR <double> p_l  (dl_size, 0);

		SEQ_VECTOR <double> s_l  (dl_size, 0);
		SEQ_VECTOR <double> q_l  (dl_size, 0);
		SEQ_VECTOR <double> m_l  (dl_size, 0);
		SEQ_VECTOR <double> n_l  (dl_size, 0);

		SEQ_VECTOR <double> u_l  (dl_size, 0);

		double gama_l  = 0;
		double gama_lp = 0;

		double delta_l = 0;
		double beta_l  = 0;

		double alpha_l = 0;
		double alpha_lp;

		double norm_l;
		double tol;


		// *** CG iteration loop ******************************************************************
		// *** - input to CG - vec_b - right hand side in primal

		// *** convert right hand side to dual
//		cluster.domains[0].multKplusLocal(in_right_hand_side_primal[0],  cluster.x_prim_cluster1[0]);
//		cluster.domains[0].B1_comp.MatVec(cluster.x_prim_cluster1[0], cluster.compressed_tmp, 'N', 0, 0, 0.0);
//		for (int d = 1; d < cluster.domains.size(); d++) {
//			cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[0]);
//			cluster.domains[d].B1_comp.MatVec(cluster.x_prim_cluster1[0], cluster.compressed_tmp, 'N', 0, 0, 1.0);
//		}
//		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);

		//// *** convert right hand side to dual
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (int d = 0; d < cluster.domains.size(); d++) {
			// *** convert right hand side to dual
			cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[d]);

			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);


		int iter = 0;
		timing.totalTime.Reset();

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);

		cilk_for (int i = 0; i < r_l.size(); i++) {	// r = b - Ax;
			r_l[i] = b_l[i] - Ax_l[i];
			wp_l[i] = 0.0;
			yp_l[i] = 0.0;
		}

		if (USE_PREC == 1) {
			apply_prec_comp_dom_B(timeEvalPrec, cluster, r_l, u_l);
		} else {
			cilk_for (int i = 0; i < r_l.size(); i++)
				u_l = r_l;
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, w_l);

		tol = epsilon * parallel_norm_compressed(cluster, b_l);

		for (iter = 0; iter < CG_max_iter; iter++) {

			timing.totalTime.AddStart(omp_get_wtime());

			alpha_lp = alpha_l;
			gama_lp  = gama_l;

			ddot_time.AddStart(omp_get_wtime());
			MPI_Request mpi_req;
			MPI_Status mpi_stat;

			SEQ_VECTOR <double> reduction_tmp (2,0);
			SEQ_VECTOR <double> send_buf (2,0);
			parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, &mpi_req, reduction_tmp, send_buf);

			ddot_time.AddEnd(omp_get_wtime());

			if (USE_PREC == 1) {

				prec_time.AddStart(omp_get_wtime());
				apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, m_l);
				prec_time.AddEnd(omp_get_wtime());

			} else {

				cilk_for (int i = 0; i < m_l.size(); i++)
					m_l[i] = w_l[i];

			}

			//------------------------------------------
			appA_time.AddStart(omp_get_wtime());
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, n_l);
			appA_time.AddEnd(omp_get_wtime());
			//------------------------------------------


#ifndef WIN32
#ifdef USE_MPI_3
			MPI_Wait(&mpi_req, &mpi_stat);
#endif
#endif
			gama_l  = reduction_tmp[0];
			delta_l = reduction_tmp[1];

			//------------------------------------------
			vec_time.AddStart(omp_get_wtime());
			if (iter == 0) {
				beta_l  = 0;
				alpha_l = gama_l / delta_l;
			} else {
				beta_l  = gama_l / gama_lp;
				alpha_l = gama_l / (delta_l - beta_l * gama_l / alpha_lp);
			}

			cilk_for (int i = 0; i < r_l.size(); i++) {
				z_l[i] = n_l[i] + beta_l  * z_l[i];
				q_l[i] = m_l[i] + beta_l  * q_l[i];
				s_l[i] = w_l[i] + beta_l  * s_l[i];
				p_l[i] = u_l[i] + beta_l  * p_l[i];
				x_l[i] = x_l[i] + alpha_l * p_l[i];
				r_l[i] = r_l[i] - alpha_l * s_l[i];
				u_l[i] = u_l[i] - alpha_l * q_l[i];
				w_l[i] = w_l[i] - alpha_l * z_l[i];
			}
			vec_time.AddEnd(omp_get_wtime());


			timing.totalTime.AddEnd(omp_get_wtime());
			//timing.totalTime.PrintLastStatMPI(0.0);

			norm_time.AddStart(omp_get_wtime());
			//norm_l = parallel_norm_compressed(cluster, u_l);
			norm_l = parallel_norm_compressed(cluster, r_l);
			norm_time.AddEnd(omp_get_wtime());

			if (mpi_rank == mpi_root) {
				//printf (       "Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
				printf ( "Iter MPI %d - norm %1.30f - tol %1.30f \r" , iter+1, norm_l, tol);
				//if (log_active == 1)
				//fprintf(stream,"Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
			}

			if (norm_l < tol)
				break;

		} // end iter loop

		// reconstruction of u

		//cilk_
		for (int d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1t_comp_dom.cols, 0.0 );

			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_l[ cluster.domains[d].lambda_map_sub_local[i]];

			cluster.domains[d].B1t_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');

			for(int i = 0; i < in_right_hand_side_primal[d].size(); i++)
				out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];

			cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);

		}


//		for (int d = 0; d < cluster.domains.size(); d++) {
//			cluster.domains[d].B1t_comp.MatVec(x_l, cluster.x_prim_cluster1[d], 'N');
//
//			cilk_for(int i = 0; i < in_right_hand_side_primal[d].size(); i++)
//				out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];
//
//			cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);
//		}
		// *** output of CG - vec_u_n in primal

}




// *** Dynamic Elasticity Solver - uses non-singular CG Solvers *******************************************
void IterSolver::Solve_Dynamic ( Cluster & cluster, string result_file, SEQ_VECTOR < SEQ_VECTOR < SEQ_VECTOR <double > > > & prim_solution_out)
{
	SEQ_VECTOR < SEQ_VECTOR <double> > vec_u     (cluster.domains.size());
	SEQ_VECTOR < SEQ_VECTOR <double> > vec_v     (cluster.domains.size());
	SEQ_VECTOR < SEQ_VECTOR <double> > vec_w     (cluster.domains.size());

    SEQ_VECTOR < SEQ_VECTOR <double> > vec_u_n   (cluster.domains.size());
	SEQ_VECTOR < SEQ_VECTOR <double> > vec_v_n   (cluster.domains.size());
	SEQ_VECTOR < SEQ_VECTOR <double> > vec_w_n   (cluster.domains.size());

	SEQ_VECTOR < SEQ_VECTOR <double> > vec_b     (cluster.domains.size());
	SEQ_VECTOR < SEQ_VECTOR <double> > vec_t_tmp (cluster.domains.size());

	for (int d = 0; d < cluster.domains.size(); d++) {
		vec_u[d]    .resize(cluster.domains[d].domain_prim_size, 0);
		vec_v[d]    .resize(cluster.domains[d].domain_prim_size, 0);
		vec_w[d]    .resize(cluster.domains[d].domain_prim_size, 0);

		vec_u_n[d]  .resize(cluster.domains[d].domain_prim_size, 0);
		vec_v_n[d]  .resize(cluster.domains[d].domain_prim_size, 0);
		vec_w_n[d]  .resize(cluster.domains[d].domain_prim_size, 0);

		vec_b[d]    .resize(cluster.domains[d].domain_prim_size, 0);
		vec_t_tmp[d].resize(cluster.domains[d].domain_prim_size, 0);
	}

	// *** Set up the initial acceleration ***********************
	for (int d = 0; d < cluster.domains.size(); d++) {
		for (int i = 2; i < vec_w[d].size(); i=i+3) {
			vec_w[d][i] = 1.0;
		}
	}
	// *** END - Set up the initial accel. ***********************


	const_beta   = cluster.dynamic_beta;
	const_deltat = cluster.dynamic_timestep;
	const_gama   = cluster.dynamic_gama;


	SEQ_VECTOR <double> const_a (8,0);
	const_a[0] = 1.0 / (const_beta * const_deltat * const_deltat);
	const_a[1] = const_gama / (const_beta * const_deltat);
	const_a[2] = 1.0 / (const_beta * const_deltat);
	const_a[3] = (1.0 / (2 * const_beta)) - 1.0;
	const_a[4] = (const_gama / const_beta) - 1.0;
	const_a[5] = const_deltat * ((const_gama / (2.0 * const_beta)) - 1.0);
	const_a[6] = const_deltat * (1.0 - const_gama);
	const_a[7] = const_deltat * const_gama;

	for (int time = 0; time < 100; time++) {
	//for (int time = 0; time < NumberOfTimeIterations; time++) {

		// *** calculate the right hand side in primal ********************************************
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			for(int i = 0; i < vec_u[d].size(); i++) {
				vec_t_tmp[d][i] = const_a[0] * vec_u[d][i] + const_a[2] * vec_v[d][i] + const_a[3] * vec_w[d][i];
			}
			cluster.domains[d].M.MatVec(vec_t_tmp[d], vec_b[d],'N');
		}

		// *** Run the CG solver **************************************************************

		if ( USE_PIPECG == 1 ) {
			//Solve_PipeCG_nonsingular( cluster, vec_b, vec_u_n);
		} else {
			Solve_RegCG_nonsingular ( cluster, vec_b, vec_u_n);
		}

		// *** END - Run the CG solver ********************************************************

		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			for(int i = 0; i < vec_u[d].size(); i++) {
				vec_w_n[d][i] = (const_a[0] * (vec_u_n[d][i] - vec_u[d][i])) - (const_a[2] * vec_v[d][i]) - (const_a[3] * vec_w  [d][i]);
				vec_v_n[d][i] = vec_v[d][i]                  + (const_a[6] * vec_w[d][i])                 + (const_a[7] * vec_w_n[d][i]);
			//}

			//for(int i = 0; i < vec_u[d].size(); i++) {
				vec_u[d][i] = vec_u_n[d][i];
				vec_v[d][i] = vec_v_n[d][i];
				vec_w[d][i] = vec_w_n[d][i];
			}
		}

		prim_solution_out.push_back(vec_u_n);

		// *** Save results to file *****************************************************
		if (FIND_SOLUTION == 1) {
			char number[8];
			sprintf(number, "%06d", time + 1);
			string result_file_w_time = string(result_file) + number;
			Save_to_Ensight_file (cluster, result_file_w_time, vec_u_n );
		}

		// *** XXX
		if (mpi_rank == mpi_root) {
			cout<<endl<< "Time iter " << time << "\t";
		}


		// *** XXX
		timing.totalTime.PrintStatMPI(0.0);
		timing.totalTime.Reset();

	} // *** END - time iter loop *******************************************************

	preproc_timing.PrintStatsMPI();
	timeEvalAppa  .PrintStatsMPI();

	if (USE_PREC == 1)
		timeEvalPrec.PrintStatsMPI();

}


// *** Coarse problem routines *******************************************
void IterSolver::CreateGGt( Cluster & cluster )      //, int mpi_rank, int mpi_root, int mpi_size, SparseSolver & GGt )

{

	double start = omp_get_wtime();

	double sc1 = omp_get_wtime();
	SparseMatrix G;


	//if (mpi_rank == mpi_root)
	//	G.MatAppend(cluster.G1);

	//for (int mr = 1; mr < mpi_size; mr++) {
	//	SparseMatrix Gtmp;
	//	SendMatrix(mpi_rank, mr, cluster.G1, mpi_root, Gtmp);

	//	if (mpi_rank == mpi_root) {
	//		G.MatAppend(Gtmp);
	//		Gtmp.Clear();
	//	}
	//}

	//// **** Log N MPI reduce
	int count_cv = 0;
	for (int li = 2; li <= 2*mpi_size; li = li * 2 ) {

		SparseMatrix recv_m;

		if (mpi_rank % li == 0) {
			if (li == 2)
				G.MatAppend(cluster.G1);
			if ((mpi_rank + li/2) < mpi_size) {
				SendMatrix(mpi_rank, mpi_rank + li/2, cluster.G1, mpi_rank,        recv_m);
				G.MatAppend(recv_m);
			} else {
				SendMatrix(mpi_rank, mpi_size + 1, cluster.G1, mpi_size + 1,        recv_m);
			}
		} else {

			if ((mpi_rank + li/2) % li == 0)
			{
				if (li == 2)
					SendMatrix(mpi_rank, mpi_rank       , cluster.G1, mpi_rank - li/2, recv_m);
				else
					SendMatrix(mpi_rank, mpi_rank       , G         , mpi_rank - li/2, recv_m);
			} else {
				SendMatrix(mpi_rank, mpi_rank+1, cluster.G1, mpi_rank+1,recv_m);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		count_cv += mpi_size/li;
		if (mpi_rank == 0)
			printf(" Collecting matrices G : %d of %d \r", count_cv, mpi_size);

	}

	//SparseMatrix Gtt;
	if (mpi_rank != mpi_root)
		G.Clear();
	//else {
	//Gtt = G;
	//G.Clear();
	//}
	// ****

	double ec1 = omp_get_wtime();


	double sp1 = omp_get_wtime();
	if (mpi_rank == mpi_root) {
		MKL_Set_Num_Threads(16);
		// Create Gt and later GGt matrices and remove all elements under main diagonal of the GGt
		SparseMatrix Gt;

		double t1 = omp_get_wtime();
		G.MatTranspose(Gt);
		cout << "Gtranspose = " << omp_get_wtime() - t1 << endl;

		t1 = omp_get_wtime();
		SparseMatrix GGt_Mat;
		GGt_Mat.MatMat(G, 'N', Gt);
		cout << "G x Gt = " << omp_get_wtime() - t1 << endl;

		t1 = omp_get_wtime();
		Gt.Clear();
		G.Clear();
		cout << "G and Gt clear = " << omp_get_wtime() - t1 << endl;

		SpyText(GGt_Mat);

		t1 = omp_get_wtime();
		GGt_Mat.RemoveLower();
		cout << "GGt remove lower = " << omp_get_wtime() - t1 << endl;

		t1 = omp_get_wtime();
		// Create Sparse Direct solver for GGt
		GGt.msglvl = 1;

		t1 = omp_get_wtime();
		GGt.ImportMatrix(GGt_Mat);
		cout << "ImportMatrix = " << omp_get_wtime() - t1 << endl;


		t1 = omp_get_wtime();
		GGt_Mat.Clear();


		t1 = omp_get_wtime();
		GGt.Factorization();
		cout << "Factorization = " << omp_get_wtime() - t1 << endl;


		t1 = omp_get_wtime();
		GGt.msglvl = 0;

		MKL_Set_Num_Threads(1);
	}

	double ep1 = omp_get_wtime();


	if (mpi_rank == mpi_root)
		GGtsize = GGt.cols;

	MPI_Bcast( & GGtsize, 1, MPI_INT, 0, MPI_COMM_WORLD);

#if TIME_MEAS >= 1
	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double end = omp_get_wtime();
		cout <<	"CG Loop - Create GGt  - collect all matrices   - Runtime = " << ec1 - sc1 << " s \n";
		cout <<	"CG Loop - Create GGt  - GGt fact. processing   - Runtime = " << ep1 - sp1 << " s \n";
		cout <<	"CG Loop - Create GGt  - total = proc + comm    - Runtime = " << end - start << " s \n\n";
	}
#endif

}

void IterSolver::CreateGGt_inv( Cluster & cluster )  //, int mpi_rank, int mpi_root, int mpi_size, SparseSolver & GGt )
{
	// temp variables
	SparseMatrix G;
	SparseMatrix Gt;
	SparseMatrix GGt_Mat_tmp;
	SparseSolver GGt_tmp;

	if (mpi_rank == 0) cout << endl;

	//// **** Log N MPI redukce
	TimeEvent collectG1_time("Collect G1 matrices to master");
	collectG1_time.AddStart(omp_get_wtime());

	int count_cv = 0;
	for (int li = 2; li <= 2*mpi_size; li = li * 2 ) {

		SparseMatrix recv_m;

		if (mpi_rank % li == 0) {
			if (li == 2)
				G.MatAppend(cluster.G1);
			if ((mpi_rank + li/2) < mpi_size) {
				SendMatrix(mpi_rank, mpi_rank + li/2, cluster.G1, mpi_rank,        recv_m);
				G.MatAppend(recv_m);
			} else {
				SendMatrix(mpi_rank, mpi_size + 1, cluster.G1, mpi_size + 1,        recv_m);
			}
		} else {

			if ((mpi_rank + li/2) % li == 0)
			{
				if (li == 2)
					SendMatrix(mpi_rank, mpi_rank       , cluster.G1, mpi_rank - li/2, recv_m);
				else
					SendMatrix(mpi_rank, mpi_rank       , G         , mpi_rank - li/2, recv_m);
			} else {
				SendMatrix(mpi_rank, mpi_rank+1, cluster.G1, mpi_rank+1,recv_m);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		count_cv += mpi_size/li;
		if (mpi_rank == 0)
			printf(" Collecting matrices G : %d of %d \r", count_cv, mpi_size);

	}
	collectG1_time.AddEnd(omp_get_wtime());
	collectG1_time.PrintStatMPI(0.0);
	preproc_timing.AddEvent(collectG1_time);

	if (mpi_rank != mpi_root)
		G.Clear();

	TimeEvent master_processing_time("Master proc.- Gtransp., GxGt mult., .");
	master_processing_time.AddStart(omp_get_wtime());

	MKL_Set_Num_Threads(24);

	if (mpi_rank == mpi_root) {
		// Create Gt and later GGt matrices and remove all elements under main diagonal of the GGt

		double t1 = omp_get_wtime();
		G.MatTranspose(Gt);
		cout << "Gtranspose = " << omp_get_wtime() - t1 << endl;

		t1 = omp_get_wtime();
		GGt_Mat_tmp.MatMat(G, 'N', Gt);
		cout << "G x Gt = " << omp_get_wtime() - t1 << endl;

		t1 = omp_get_wtime();
		Gt.Clear();
		G.Clear();
		cout << "G and Gt clear = " << omp_get_wtime() - t1 << endl;

		t1 = omp_get_wtime();
		GGt_Mat_tmp.RemoveLower();
		cout << "GGt remove lower = " << omp_get_wtime() - t1 << endl;


		SpyText(GGt_Mat_tmp);


	}
	master_processing_time.AddEnd(omp_get_wtime());
	master_processing_time.PrintStatMPI(0.0);
	preproc_timing.AddEvent(master_processing_time);

	 TimeEvent GGt_bcast_time("Time to broadcast GGt from master all");
	 GGt_bcast_time.AddStart(omp_get_wtime());
	BcastMatrix(mpi_rank, mpi_root, mpi_root, GGt_Mat_tmp);
	 GGt_bcast_time.AddEnd(omp_get_wtime());
	 GGt_bcast_time.PrintStatMPI(0.0);
	 preproc_timing.AddEvent(GGt_bcast_time);


	// Create Sparse Direct solver for GGt
	if (mpi_rank == mpi_root)
		GGt_tmp.msglvl = 1;

	 TimeEvent importGGt_time("Time to import GGt matrix into solver");
	 importGGt_time.AddStart(omp_get_wtime());
	GGt_tmp.ImportMatrix(GGt_Mat_tmp);
	 importGGt_time.AddEnd(omp_get_wtime());
	 importGGt_time.PrintStatMPI(0.0);
	 preproc_timing.AddEvent(importGGt_time);

	GGt_Mat_tmp.Clear();

	 TimeEvent GGtFactor_time("GGT Factorization time");
	 GGtFactor_time.AddStart(omp_get_wtime());
	GGt_tmp.Factorization();
	 GGtFactor_time.AddEnd(omp_get_wtime());
	 GGtFactor_time.PrintStatMPI(0.0);
	 preproc_timing.AddEvent(GGtFactor_time);

	TimeEvent GGT_rhs_time("Time to create RHS for get GGTINV");
	GGT_rhs_time.AddStart(omp_get_wtime());
	SEQ_VECTOR <double> rhs   (cluster.G1.rows * GGt_tmp.rows, 0);
	cluster.GGtinvV.resize(cluster.G1.rows * GGt_tmp.rows, 0);

	for (int i = 0; i < cluster.G1.rows; i++) {
		int index = (GGt_tmp.rows * i) + (cluster.G1.rows * mpi_rank) + i;
		rhs[index] = 1;
	}
	 GGT_rhs_time.AddEnd(omp_get_wtime());
	 GGT_rhs_time.PrintStatMPI(0.0);
	 preproc_timing.AddEvent(GGT_rhs_time);

	 TimeEvent GGt_solve_time("Running solve to get stripe(s) of GGtINV");
	 GGt_solve_time.AddStart(omp_get_wtime());

	GGt_tmp.Solve(rhs, cluster.GGtinvV, cluster.G1.rows);

	cluster.GGtinvM.dense_values = cluster.GGtinvV;
	cluster.GGtinvM.cols = cluster.G1.rows;
	cluster.GGtinvM.rows = GGt_tmp.rows;

	GGtsize  = GGt_tmp.cols;

	GGt.cols = GGt_tmp.cols;
	GGt.rows = GGt_tmp.rows;
	GGt.nnz  = GGt_tmp.nnz;

	GGt_tmp.msglvl = 0;
	GGt_tmp.Clear();

	 GGt_solve_time.AddEnd(omp_get_wtime());
	 GGt_solve_time.PrintStatMPI(0.0);
	 preproc_timing.AddEvent(GGt_solve_time);

	 MKL_Set_Num_Threads(1);

}


void IterSolver::CreateGGt_inv_dist( Cluster & cluster )  //, int mpi_rank, int mpi_root, int mpi_size, SparseSolver & GGt )
{

	// temp variables
	vector < SparseMatrix > G_neighs   ( cluster.my_neighs.size() );
	vector < SparseMatrix > GGt_neighs ( cluster.my_neighs.size() );
	SparseMatrix Gt_l;
	SparseMatrix GGt_l;
	SparseMatrix GGt_Mat_tmp;
	SparseSolver GGt_tmp;

	 TimeEvent SaRGlocal("Send a Receive local G1 matrices to neighs. "); SaRGlocal.AddStart(omp_get_wtime());
	for (int neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ )
		SendMatrix(cluster.G1, cluster.my_neighs[neigh_i]);

	for (int neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ )
		RecvMatrix(G_neighs[neigh_i], cluster.my_neighs[neigh_i]);
	 SaRGlocal.AddEnd(omp_get_wtime()); SaRGlocal.PrintStatMPI(0.0); preproc_timing.AddEvent(SaRGlocal);


	 TimeEvent Gt_l_trans("Local G1 matrix transpose to create Gt "); Gt_l_trans.AddStart(omp_get_wtime());
	if (cluster.USE_HFETI == 0)
		cluster.G1.MatTranspose(Gt_l);
	 Gt_l_trans.AddEnd(omp_get_wtime()); Gt_l_trans.PrintStatMPI(0.0); preproc_timing.AddEvent(Gt_l_trans);

	 TimeEvent GxGtMatMat("Local G x Gt MatMat "); GxGtMatMat.AddStart(omp_get_wtime());

	if (cluster.USE_HFETI == 0)
		GGt_l.MatMat(cluster.G1, 'N', Gt_l);
	else
		GGt_l.MatMatT(cluster.G1, cluster.G1);

	 GxGtMatMat.AddEnd(omp_get_wtime()); GxGtMatMat.PrintStatMPI(0.0); preproc_timing.AddEvent(GxGtMatMat);
	 //GxGtMatMat.PrintLastStatMPI_PerNode(0.0);


	for (int i = 0; i < GGt_l.CSR_J_col_indices.size(); i++)
		GGt_l.CSR_J_col_indices[i] += mpi_rank * cluster.G1.rows;

	GGt_l.cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;

	 TimeEvent GGTNeighTime("G1t_local x G1_neigh MatMat(N-times) "); GGTNeighTime.AddStart(omp_get_wtime());
	cilk_for (int neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {

		if (cluster.USE_HFETI == 0)
			GGt_neighs[neigh_i].MatMat(G_neighs[neigh_i], 'N', Gt_l);
		else
			GGt_neighs[neigh_i].MatMatT(G_neighs[neigh_i], cluster.G1);

		GGt_neighs[neigh_i].MatTranspose();

		int inc = cluster.G1.rows * cluster.my_neighs[neigh_i];
		for (int i = 0; i < GGt_neighs[neigh_i].CSR_J_col_indices.size(); i++)
			GGt_neighs[neigh_i].CSR_J_col_indices[i] += inc;

		GGt_neighs[neigh_i].cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;
		G_neighs[neigh_i].Clear();
	}
	 GGTNeighTime.AddEnd(omp_get_wtime()); GGTNeighTime.PrintStatMPI(0.0); preproc_timing.AddEvent(GGTNeighTime);
	 //GGTNeighTime.PrintLastStatMPI_PerNode(0.0);

	 TimeEvent GGtLocAsm("Assembling row of GGt per node - MatAddInPlace "); GGtLocAsm.AddStart(omp_get_wtime());
	for (int neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		GGt_l.MatAddInPlace(GGt_neighs[neigh_i], 'N', 1.0);
		GGt_neighs[neigh_i].Clear();
	}
	 GGtLocAsm.AddEnd(omp_get_wtime()); GGtLocAsm.PrintStatMPI(0.0); preproc_timing.AddEvent(GGtLocAsm);


	TimeEvent collectGGt_time("Collect GGt pieces to master"); 	collectGGt_time.AddStart(omp_get_wtime());
	int count_cv_l = 0;
	for (int li = 2; li <= 2*mpi_size; li = li * 2 ) {

		SparseMatrix recv_m_l;

		if (mpi_rank % li == 0) {
			if (li == 2)
				GGt_Mat_tmp.MatAppend(GGt_l);
			if ((mpi_rank + li/2) < mpi_size) {
				SendMatrix(mpi_rank, mpi_rank + li/2, GGt_l, mpi_rank,     recv_m_l);
				GGt_Mat_tmp.MatAppend(recv_m_l);
			} else {
				SendMatrix(mpi_rank, mpi_size + 1   , GGt_l, mpi_size + 1, recv_m_l);
			}
		} else {

			if ((mpi_rank + li/2) % li == 0)
			{
				if (li == 2)
					SendMatrix(mpi_rank, mpi_rank       , GGt_l      , mpi_rank - li/2, recv_m_l);
				else
					SendMatrix(mpi_rank, mpi_rank       , GGt_Mat_tmp, mpi_rank - li/2, recv_m_l);
			} else {
				SendMatrix(mpi_rank, mpi_rank+1, GGt_l, mpi_rank+1,recv_m_l);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		GGt_l.Clear();

		count_cv_l += mpi_size/li;
		if (mpi_rank == 0)
			printf(" Collecting matrices G : %d of %d \r", count_cv_l, mpi_size);

	}
	collectGGt_time.AddEnd(omp_get_wtime()); collectGGt_time.PrintStatMPI(0.0); preproc_timing.AddEvent(collectGGt_time);

	if (mpi_rank == 0) cout << endl;

	if (mpi_rank == 0)  {
		GGt_Mat_tmp.RemoveLower();
		SpyText(GGt_Mat_tmp);
	}

	MKL_Set_Num_Threads(24);

	TimeEvent GGt_bcast_time("Time to broadcast GGt from master all"); GGt_bcast_time.AddStart(omp_get_wtime());
	BcastMatrix(mpi_rank, mpi_root, mpi_root, GGt_Mat_tmp);
	GGt_bcast_time.AddEnd(omp_get_wtime()); GGt_bcast_time.PrintStatMPI(0.0); preproc_timing.AddEvent(GGt_bcast_time);

	// Create Sparse Direct solver for GGt
	if (mpi_rank == mpi_root)
		GGt_tmp.msglvl = 1;

	TimeEvent importGGt_time("Time to import GGt matrix into solver"); importGGt_time.AddStart(omp_get_wtime());
	GGt_tmp.ImportMatrix(GGt_Mat_tmp);
	importGGt_time.AddEnd(omp_get_wtime()); importGGt_time.PrintStatMPI(0.0); preproc_timing.AddEvent(importGGt_time);

	GGt_Mat_tmp.Clear();

	TimeEvent GGtFactor_time("GGT Factorization time"); GGtFactor_time.AddStart(omp_get_wtime());
	GGt_tmp.SetThreaded();
	GGt_tmp.Factorization();
	GGtFactor_time.AddEnd(omp_get_wtime()); GGtFactor_time.PrintStatMPI(0.0); preproc_timing.AddEvent(GGtFactor_time);

	TimeEvent GGT_rhs_time("Time to create RHS for get GGTINV"); GGT_rhs_time.AddStart(omp_get_wtime());
	SEQ_VECTOR <double> rhs   (cluster.G1.rows * GGt_tmp.rows, 0);
	cluster.GGtinvV.resize(cluster.G1.rows * GGt_tmp.rows, 0);

	for (int i = 0; i < cluster.G1.rows; i++) {
		int index = (GGt_tmp.rows * i) + (cluster.G1.rows * mpi_rank) + i;
		rhs[index] = 1;
	}
	GGT_rhs_time.AddEnd(omp_get_wtime()); GGT_rhs_time.PrintStatMPI(0.0); preproc_timing.AddEvent(GGT_rhs_time);

	TimeEvent GGt_solve_time("Running solve to get stripe(s) of GGtINV"); GGt_solve_time.AddStart(omp_get_wtime());

	GGt_tmp.Solve(rhs, cluster.GGtinvV, cluster.G1.rows);

	cluster.GGtinvM.dense_values = cluster.GGtinvV;
	cluster.GGtinvM.cols = cluster.G1.rows;
	cluster.GGtinvM.rows = GGt_tmp.rows;

	GGtsize  = GGt_tmp.cols;

	GGt.cols = GGt_tmp.cols;
	GGt.rows = GGt_tmp.rows;
	GGt.nnz  = GGt_tmp.nnz;

	GGt_tmp.msglvl = 0;
	GGt_tmp.Clear();

	GGt_solve_time.AddEnd(omp_get_wtime()); GGt_solve_time.PrintStatMPI(0.0); preproc_timing.AddEvent(GGt_solve_time);

	MKL_Set_Num_Threads(1);

}


// *** Projector routines ************************************************
void IterSolver::Projector_l_compG (TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, int output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // int mpi_rank, SparseSolver & GGt,
{

	time_eval.totalTime.AddStart(omp_get_wtime());

	//int dual_size    = cluster.domains[0].B1.rows;
	int d_local_size = cluster.G1_comp.rows;
	int mpi_root     = 0;

	SEQ_VECTOR<double> d_local( d_local_size );
	SEQ_VECTOR<double> d_mpi  ( GGtsize );

	time_eval.timeEvents[0].AddStart(omp_get_wtime());
	if ( output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 1)
		d_local = x_in;
	else
		cluster.G1_comp.MatVec(x_in, d_local, 'N');
	time_eval.timeEvents[0].AddEnd(omp_get_wtime());



	time_eval.timeEvents[1].AddStart(omp_get_wtime());
	MPI_Gather(&d_local[0], d_local_size, MPI_DOUBLE,
		&d_mpi[0], d_local_size, MPI_DOUBLE,
		mpi_root, MPI_COMM_WORLD);
	time_eval.timeEvents[1].AddEnd(omp_get_wtime());




	time_eval.timeEvents[2].AddStart(omp_get_wtime());
	if (mpi_rank == mpi_root ) {
		GGt.Solve(d_mpi);				// t1 = Uc\(Lc\d);
	}
	time_eval.timeEvents[2].AddEnd(omp_get_wtime());




	time_eval.timeEvents[3].AddStart(omp_get_wtime());
	MPI_Scatter( &d_mpi[0],      d_local_size, MPI_DOUBLE,
		&d_local[0], d_local_size, MPI_DOUBLE,
		mpi_root, MPI_COMM_WORLD);
	time_eval.timeEvents[3].AddEnd(omp_get_wtime());

	if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 2) {
		// for mu calculation
		y_out = d_local;

	} else {

		time_eval.timeEvents[4].AddStart(omp_get_wtime());
		//cluster.G1t_comp.MatVec(d_local, cluster.compressed_tmp, 'N'); // SUPER POZOR
		cluster.G1_comp.MatVec(d_local, cluster.compressed_tmp, 'T');
		time_eval.timeEvents[4].AddEnd(omp_get_wtime());

		time_eval.timeEvents[5].AddStart(omp_get_wtime());
		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
		time_eval.timeEvents[5].AddEnd(omp_get_wtime());

		if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 0) {
			cilk_for (int i = 0; i < x_in.size(); i++)
				y_out[i] = x_in[i] - y_out[i];
		}

	}

	time_eval.totalTime.AddEnd(omp_get_wtime());

}

void IterSolver::Projector_l_inv_compG (TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, int output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // int mpi_rank, SparseSolver & GGt,
{

	time_eval.totalTime.AddStart(omp_get_wtime());

	int d_local_size = cluster.G1_comp.rows;
	int mpi_root     = 0;

	SEQ_VECTOR<double> d_local( d_local_size );
	SEQ_VECTOR<double> d_mpi  ( GGtsize );

	time_eval.timeEvents[0].AddStart(omp_get_wtime());
	if ( output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 1)
		d_local = x_in;
	else
		cluster.G1_comp.MatVec(x_in, d_local, 'N');
	time_eval.timeEvents[0].AddEnd(omp_get_wtime());

	time_eval.timeEvents[1].AddStart(omp_get_wtime());
	MPI_Allgather(&d_local[0], d_local_size, MPI_DOUBLE,
		&d_mpi[0], d_local_size, MPI_DOUBLE,
		MPI_COMM_WORLD);
	time_eval.timeEvents[1].AddEnd(omp_get_wtime());

	time_eval.timeEvents[2].AddStart(omp_get_wtime());
	//for (int j = 0; j < cluster.G1_comp.rows; j++) {
	//	d_local[j] = 0;
	//	for (int i = 0; i < GGt.rows; i++ ) {
	//		d_local[j] += cluster.GGtinvV[j * GGt.rows + i] * d_mpi[i];			// t1 = Uc\(Lc\d);
	//	}
	//}
	cluster.GGtinvM.DenseMatVec(d_mpi, d_local, 'T');
	time_eval.timeEvents[2].AddEnd(omp_get_wtime());

	time_eval.timeEvents[3].AddStart(omp_get_wtime());
	//MPI_Scatter( &d_mpi[0],      d_local_size, MPI_DOUBLE,
	//	&d_local[0], d_local_size, MPI_DOUBLE,
	//	mpi_root, MPI_COMM_WORLD);
	time_eval.timeEvents[3].AddEnd(omp_get_wtime());

	if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 2) {
		y_out = d_local; // for RBM amplitudes calculation
	} else {

		time_eval.timeEvents[4].AddStart(omp_get_wtime());
		//cluster.G1t_comp.MatVec(d_local, cluster.compressed_tmp, 'N'); SUPER POZOR
		cluster.G1_comp.MatVec(d_local, cluster.compressed_tmp, 'T');
		time_eval.timeEvents[4].AddEnd(omp_get_wtime());

		time_eval.timeEvents[5].AddStart(omp_get_wtime());
		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
		time_eval.timeEvents[5].AddEnd(omp_get_wtime());

		if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 0) {
			cilk_for (int i = 0; i < y_out.size(); i++)
				y_out[i] = x_in[i] - y_out[i];
		}

	}

	time_eval.totalTime.AddEnd(omp_get_wtime());

}

// *** Action of K+ routines *********************************************
void IterSolver::apply_A_l_compB( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {

	time_eval.totalTime.AddStart(omp_get_wtime());

	time_eval.timeEvents[0].AddStart(omp_get_wtime());
	if (cluster.USE_KINV == 0) {
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			cluster.domains[d].B1t_comp.MatVec (x_in, cluster.x_prim_cluster1[d], 'N'); //cluster.domains[d].B1_comp.MatVec(x_in, cluster.x_prim_cluster1[d], 'T'); //cluster.domains[d].B1_comp.MatVecCOO(x_in, cluster.x_prim_cluster1[d], 'T');
		}
	}
	time_eval.timeEvents[0].AddEnd(omp_get_wtime());


	time_eval.timeEvents[1].AddStart(omp_get_wtime());
	if (cluster.USE_HFETI == 0) {
		if (cluster.USE_KINV == 0) {

			cilk_for (int d = 0; d < cluster.domains.size(); d++)
				cluster.domains[d].multKplusLocal(cluster.x_prim_cluster1[d]);

		} else {

			cilk_for (int d = 0; d < cluster.domains.size(); d++) {
				if (d == 0)
					cluster.domains[0].B1Kplus.DenseMatVec(x_in, cluster.compressed_tmp);
				else
					cluster.domains[d].B1Kplus.DenseMatVec(x_in, cluster.domains[d].compressed_tmp);
			}

			for ( int d = 1; d < cluster.domains.size(); d++ )
				for ( int i = 0; i < cluster.domains[d].compressed_tmp.size(); i++ )
					cluster.compressed_tmp[i] = cluster.compressed_tmp[i] + cluster.domains[d].compressed_tmp[i];

		}
	} else {
		if (cluster.USE_KINV == 0){
			cluster.multKplusGlobal_l(cluster.x_prim_cluster1);
		} else {
			cout << "USE of K_INV is not supported for HTFETI";
			exit(0);
		}
	}
	time_eval.timeEvents[1].AddEnd(omp_get_wtime());


	time_eval.timeEvents[2].AddStart(omp_get_wtime());
	if (cluster.USE_KINV == 0) {
		cluster.B1_comp_MatVecSum( cluster.x_prim_cluster1, cluster.compressed_tmp, 'N' );
	}
	time_eval.timeEvents[2].AddEnd(omp_get_wtime());


	time_eval.timeEvents[3].AddStart(omp_get_wtime());
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[3].AddEnd(omp_get_wtime());


	time_eval.totalTime.AddEnd(omp_get_wtime());

}

void IterSolver::apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {

	if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {
		time_eval.timeEvents[0].AddStart(omp_get_wtime());
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
				x_in_tmp[i]                            = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
#ifdef CUDA
				cluster.domains[d].cuda_pinned_buff[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
#endif
				//cluster.domains[d].cuda_pinned_buff_fl[i] = (float) x_in[ cluster.domains[d].lambda_map_sub_local[i]];
			}

#ifdef CUDA
			cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start( cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff,'N',0 );
			//cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start_fl( cluster.domains[d].cuda_pinned_buff_fl, cluster.domains[d].cuda_pinned_buff_fl,'N',0 );
#else
			cluster.domains[d].B1Kplus.DenseMatVec (x_in_tmp, cluster.domains[d].compressed_tmp);
#endif

			cluster.domains[d].B1t_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
		}
		time_eval.timeEvents[0].AddEnd(omp_get_wtime());

		time_eval.timeEvents[1].AddStart(omp_get_wtime());
		cluster.multKplusGlobal_Kinv( cluster.x_prim_cluster1 );
		time_eval.timeEvents[1].AddEnd(omp_get_wtime());

		time_eval.timeEvents[2].AddStart(omp_get_wtime());
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);

#ifdef CUDA
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_sync ( );
		}

		for (int d = 0; d < cluster.domains.size(); d++) {
			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].cuda_pinned_buff[i];
				//cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += (double)cluster.domains[d].cuda_pinned_buff_fl[i];
			}
		}

#else
		for (int d = 0; d < cluster.domains.size(); d++) {
			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
		}
#endif

		//std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (int d = 0; d < cluster.domains.size(); d++) {
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		time_eval.timeEvents[2].AddEnd(omp_get_wtime());

	}
	time_eval.totalTime.AddStart(omp_get_wtime());


	if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 0) {
		time_eval.timeEvents[0].AddStart(omp_get_wtime());
		//// POZOR - jen pro porovnani vykonu CPU a GPU
		//cilk_for (int d = 0; d < cluster.domains.size(); d++) {
		//	SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
		//	for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
		//		x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
		//	cluster.domains[d].B1Kplus.DenseMatVec(x_in_tmp, cluster.domains[d].compressed_tmp);
		//}
		//// POZOR - jen pro porovnani vykonu CPU a GPU
		time_eval.timeEvents[0].AddEnd(omp_get_wtime());

		time_eval.timeEvents[1].AddStart(omp_get_wtime());
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
#ifdef CUDA
			cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy(x_in_tmp, cluster.domains[d].compressed_tmp,'N',0);
#else
			cluster.domains[d].B1Kplus.DenseMatVec(x_in_tmp, cluster.domains[d].compressed_tmp);
#endif
		}
		time_eval.timeEvents[1].AddEnd(omp_get_wtime());

		time_eval.timeEvents[2].AddStart(omp_get_wtime());
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		for (int d = 0; d < cluster.domains.size(); d++) {
			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
		}
		time_eval.timeEvents[2].AddEnd(omp_get_wtime());
	}



	if (cluster.USE_KINV == 0) {
		time_eval.timeEvents[0].AddStart(omp_get_wtime());
		cilk_for (int d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1t_comp_dom.cols, 0.0 );
			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
			cluster.domains[d].B1t_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
			//cluster.x_prim_cluster2[d] = cluster.x_prim_cluster1[d]; // POZOR zbytecne kopirovani // prim norm
		}
		time_eval.timeEvents[0].AddEnd(omp_get_wtime());

		time_eval.timeEvents[1].AddStart(omp_get_wtime());
		if (cluster.USE_HFETI == 0) {
			cilk_for (int d = 0; d < cluster.domains.size(); d++)
				cluster.domains[d].multKplusLocal(cluster.x_prim_cluster1[d]);
			} else {
				cluster.multKplusGlobal_l(cluster.x_prim_cluster1);
		}
		time_eval.timeEvents[1].AddEnd(omp_get_wtime());


		time_eval.timeEvents[2].AddStart(omp_get_wtime());
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (int d = 0; d < cluster.domains.size(); d++) {
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		time_eval.timeEvents[2].AddEnd(omp_get_wtime());

	}

	time_eval.timeEvents[3].AddStart(omp_get_wtime());
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[3].AddEnd(omp_get_wtime());

	time_eval.totalTime.AddEnd(omp_get_wtime());

}

void IterSolver::apply_prec_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ) {

	time_eval.totalTime.AddStart(omp_get_wtime());

	time_eval.timeEvents[0].AddStart(omp_get_wtime());

	cilk_for (int d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1t_comp_dom.cols, 0.0 );
		for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling
		cluster.domains[d].B1t_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');

		cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');

		//cluster.x_prim_cluster2[d] = cluster.x_prim_cluster1[d];

	}

	std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
	SEQ_VECTOR < double > y_out_tmp;
	for (int d = 0; d < cluster.domains.size(); d++) {
		y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
		cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

		for (int i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling
	}
	time_eval.timeEvents[0].AddEnd(omp_get_wtime());


	time_eval.timeEvents[1].AddStart(omp_get_wtime());
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[1].AddEnd(omp_get_wtime());


	time_eval.totalTime.AddEnd(omp_get_wtime());

}


void IterSolver::apply_prec_compB( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ) {

	time_eval.totalTime.AddStart(omp_get_wtime());

	time_eval.timeEvents[0].AddStart(omp_get_wtime());
	//vector <double>  y_tmp (x_in.size(), 0);
	//for (int d = 0; d < cluster.domains.size(); d++) {
	//	vector <double > tmp (cluster.domains[d].domain_prim_size);
	//	vector <double > tmp2 (cluster.domains[d].domain_prim_size);
	//	cluster.domains[d].B1_comp.MatVec(x_in,tmp,'T');
	//	cluster.domains[d].Prec.MatVec(tmp, tmp2,'N');
	//	cluster.domains[d].B1_comp.MatVec(tmp2, y_tmp, 'N', 0, 0, 1.0); // will add (summation per elements) all partial results into y_out
	//}
	//
	//All_Reduce_lambdas_compB(cluster, y_tmp, y_out);

	cilk_for (int d = 0; d < cluster.domains.size(); d++) {
		//cluster.domains[d].B1_comp.MatVec(x_in,cluster.x_prim_cluster1[d],'T');
		cluster.domains[d].B1t_comp.MatVec(x_in,cluster.x_prim_cluster1[d],'N');
		cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
	}

	cluster.B1_comp_MatVecSum( cluster.x_prim_cluster2, cluster.compressed_tmp, 'N' );

	//cluster.domains[0].B1_comp.MatVec    (cluster.x_prim_cluster2[0], cluster.compressed_tmp, 'N', 0, 0, 0.0); // first vector overwrites temp vector
	//for (int d = 1; d < cluster.domains.size(); d++) { // reduction
	//	cluster.domains[d].B1_comp.MatVec(cluster.x_prim_cluster2[d], cluster.compressed_tmp, 'N', 0, 0, 1.0); // will add (summation per elements) all partial results into y_out
	//}

	time_eval.timeEvents[0].AddEnd(omp_get_wtime());


	time_eval.timeEvents[1].AddStart(omp_get_wtime());
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[1].AddEnd(omp_get_wtime());


	time_eval.totalTime.AddEnd(omp_get_wtime());

}

// *** END - Iteration solver class *************************************
// **********************************************************************



// **********************************************************************
// *** Communication layer **********************************************
void   SendMatrix  ( int rank, int source_rank, SparseMatrix & A_in, int dest_rank, SparseMatrix & B_out) {

	int param_tag = 1;
	int I_row_tag = 2;
	int J_col_tag = 3;
	int V_val_tag = 4;

	MPI_Status status;
	MPI_Request request;

	if (rank == source_rank) {
		int send_par_buf[4];
		send_par_buf[0] = A_in.cols;
		send_par_buf[1] = A_in.rows;
		send_par_buf[2] = A_in.nnz;
		send_par_buf[3] = A_in.type;

#ifdef XE6
		MPI_Send(send_par_buf, 4, MPI_INT, dest_rank, param_tag, MPI_COMM_WORLD);
		MPI_Send(&A_in.CSR_I_row_indices[0], A_in.rows + 1, MPI_INT, dest_rank, I_row_tag, MPI_COMM_WORLD );
		MPI_Send(&A_in.CSR_J_col_indices[0], A_in.nnz,      MPI_INT, dest_rank, J_col_tag, MPI_COMM_WORLD );
		MPI_Send(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD );
#else
		MPI_Isend(send_par_buf, 4, MPI_INT, dest_rank, param_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, MPI_INT, dest_rank, I_row_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      MPI_INT, dest_rank, J_col_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD, & request);
#endif

	}

	if (rank == dest_rank) {
		int recv_par_buf[4];
		MPI_Recv(recv_par_buf, 4, MPI_INT, source_rank, param_tag, MPI_COMM_WORLD, & status);
		B_out.cols = recv_par_buf[0];
		B_out.rows = recv_par_buf[1];
		B_out.nnz  = recv_par_buf[2];
		B_out.type = recv_par_buf[3];

		B_out.CSR_I_row_indices.resize(B_out.rows + 1);
		B_out.CSR_J_col_indices.resize(B_out.nnz);
		B_out.CSR_V_values.     resize(B_out.nnz);

		MPI_Recv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, MPI_INT,    source_rank, I_row_tag, MPI_COMM_WORLD, & status );
		MPI_Recv(&B_out.CSR_J_col_indices[0], B_out.nnz,      MPI_INT,    source_rank, J_col_tag, MPI_COMM_WORLD, & status );
		MPI_Recv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE, source_rank, V_val_tag, MPI_COMM_WORLD, & status );
	}

#ifdef WIN32
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void   SendMatrix ( SparseMatrix & A_in, int dest_rank ) {

	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);


	int param_tag = 1;
	int I_row_tag = 2;
	int J_col_tag = 3;
	int V_val_tag = 4;

	MPI_Status status;
	MPI_Request request;

	int send_par_buf[4];
	send_par_buf[0] = A_in.cols;
	send_par_buf[1] = A_in.rows;
	send_par_buf[2] = A_in.nnz;
	send_par_buf[3] = A_in.type;

//#ifdef XE6
//		MPI_Send(send_par_buf, 4, MPI_INT, dest_rank, param_tag, MPI_COMM_WORLD);
//		MPI_Send(&A_in.CSR_I_row_indices[0], A_in.rows + 1, MPI_INT, dest_rank, I_row_tag, MPI_COMM_WORLD );
//		MPI_Send(&A_in.CSR_J_col_indices[0], A_in.nnz,      MPI_INT, dest_rank, J_col_tag, MPI_COMM_WORLD );
//		MPI_Send(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD );
//#else
		MPI_Isend(send_par_buf, 4, MPI_INT, dest_rank, param_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, MPI_INT, dest_rank, I_row_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      MPI_INT, dest_rank, J_col_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD, & request);
//#endif

#ifdef WIN32
//	MPI_Barrier(MPI_COMM_WORLD);
#endif
}



void   RecvMatrix ( SparseMatrix & B_out, int source_rank) {

	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);


	int param_tag = 1;
	int I_row_tag = 2;
	int J_col_tag = 3;
	int V_val_tag = 4;

	MPI_Status status;
	MPI_Request request;

	int recv_par_buf[4];
	MPI_Recv(recv_par_buf, 4, MPI_INT, source_rank, param_tag, MPI_COMM_WORLD, & status);
	B_out.cols = recv_par_buf[0];
	B_out.rows = recv_par_buf[1];
	B_out.nnz  = recv_par_buf[2];
	B_out.type = recv_par_buf[3];

	B_out.CSR_I_row_indices.resize(B_out.rows + 1);
	B_out.CSR_J_col_indices.resize(B_out.nnz);
	B_out.CSR_V_values.     resize(B_out.nnz);

	MPI_Recv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, MPI_INT,    source_rank, I_row_tag, MPI_COMM_WORLD, & status );
	MPI_Recv(&B_out.CSR_J_col_indices[0], B_out.nnz,      MPI_INT,    source_rank, J_col_tag, MPI_COMM_WORLD, & status );
	MPI_Recv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE, source_rank, V_val_tag, MPI_COMM_WORLD, & status );


#ifdef WIN32
	//MPI_Barrier(MPI_COMM_WORLD);
#endif
}



void   BcastMatrix ( int rank, int mpi_root, int source_rank, SparseMatrix & A) {

	//int param_tag = 1;
	//int I_row_tag = 2;
	//int J_col_tag = 3;
	//int V_val_tag = 4;

	//MPI_Status status;
	//MPI_Request request;

	int send_par_buf[4];
	//int recv_par_buf[4];

	if (rank == source_rank) {
		send_par_buf[0] = A.cols; send_par_buf[1] = A.rows; send_par_buf[2] = A.nnz; send_par_buf[3] = A.type;
	}

	MPI_Bcast(send_par_buf, 4, MPI_INT, source_rank, MPI_COMM_WORLD);

	if (rank != source_rank) {
		A.cols = send_par_buf[0]; A.rows = send_par_buf[1]; A.nnz  = send_par_buf[2]; A.type = send_par_buf[3];
		A.CSR_I_row_indices.resize(A.rows + 1);
		A.CSR_J_col_indices.resize(A.nnz);
		A.CSR_V_values.     resize(A.nnz);
	}

	MPI_Bcast(&A.CSR_I_row_indices[0], A.rows + 1, MPI_INT, source_rank, MPI_COMM_WORLD);
	MPI_Bcast(&A.CSR_J_col_indices[0], A.nnz,      MPI_INT, source_rank, MPI_COMM_WORLD);
	MPI_Bcast(&A.CSR_V_values[0],      A.nnz,   MPI_DOUBLE, source_rank, MPI_COMM_WORLD);
}

void   All_Reduce_lambdas_compB( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out )
{
	for (int i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (int j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			cluster.my_comm_lambdas[i][j] = x_in[cluster.my_comm_lambdas_indices_comp[i][j]];
		}
	}


	MPI_Request * mpi_req  = new MPI_Request [cluster.my_neighs.size()];
	MPI_Status  * mpi_stat = new MPI_Status  [cluster.my_neighs.size()];

	cluster.iter_cnt_comm++;
	int tag = cluster.iter_cnt_comm;

	for (int neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		MPI_Sendrecv(
			&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,
			&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag, MPI_COMM_WORLD, &mpi_stat[neigh_i] );
	}

	//for (int neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
	//	int b_size = cluster.my_comm_lambdas[neigh_i].size();
	//	MPI_Isend(&b_size,                              1                                      , MPI_INT   , cluster.my_neighs[neigh_i], tag + 100, MPI_COMM_WORLD, &mpi_req[neigh_i] );
	//	MPI_Isend(&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,       MPI_COMM_WORLD, &mpi_req[neigh_i] );
	//
	//}

	//for (int neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
	//	int r_size = 0;
	//	MPI_Recv(&r_size                             ,                                       1, MPI_INT   , cluster.my_neighs[neigh_i], tag + 100, MPI_COMM_WORLD, &mpi_stat[neigh_i] );
	//	if (r_size != cluster.my_recv_lambdas[neigh_i].size()) cout << "Error - different buffer size " << endl;
	//	MPI_Recv(&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag      , MPI_COMM_WORLD, &mpi_stat[neigh_i] );
	//}

#ifdef XE6
	MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef WIN32
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	delete [] mpi_req;
	delete [] mpi_stat;

	y_out = x_in; // POZOR pozor
	for (int i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (int j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			y_out[cluster.my_comm_lambdas_indices_comp[i][j]] += cluster.my_recv_lambdas[i][j];
		}
	}



}

void   compress_lambda_vector  ( Cluster & cluster, SEQ_VECTOR <double> & decompressed_vec_lambda)
{
	//compress vector for CG in main loop
	for (int i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[cluster.my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(cluster.my_lamdas_indices.size());
}

void   decompress_lambda_vector( Cluster & cluster, SEQ_VECTOR <double> & compressed_vec_lambda)
{
	SEQ_VECTOR <double> decompressed_vec_lambda (cluster.domains[0].B1.rows,0);

	for (int i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[cluster.my_lamdas_indices[i]] = compressed_vec_lambda[i];

	compressed_vec_lambda = decompressed_vec_lambda;
}

double parallel_norm_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector )
{

	double wl = 0; double wg = 0;

	for (int i = 0; i < cluster.my_lamdas_indices.size(); i++)
		wl = wl + (input_vector[i] * input_vector[i] * cluster.my_lamdas_ddot_filter[i]);

	MPI_Allreduce( &wl, &wg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double norm_l = sqrt(wg);

	return norm_l;
}

double parallel_ddot_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector1, SEQ_VECTOR<double> & input_vector2 )
{
	double a1 = 0; double a1g = 0;

	for (int i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
		a1 = a1 + (input_vector1[i] * input_vector2[i] * cluster.my_lamdas_ddot_filter[i]);
	}

	MPI_Allreduce( &a1, &a1g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return a1g;
}

void   parallel_ddot_compressed_non_blocking( Cluster & cluster,
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b,
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b,
	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf)
{

	for (int i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
		send_buf[0] = send_buf[0] + (input_vector_1a[i] * input_vector_1b[i] * cluster.my_lamdas_ddot_filter[i]);
		send_buf[1] = send_buf[1] + (input_vector_2a[i] * input_vector_2b[i] * cluster.my_lamdas_ddot_filter[i]);
	}


#ifdef WIN32
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce( &send_buf[0], &output[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
#ifdef USE_MPI_3
	MPI_Iallreduce( &send_buf[0], &output[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_req);
#else
	MPI_Allreduce( &send_buf[0], &output[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif

}

// *** END - Communication layer ****************************************
// **********************************************************************
