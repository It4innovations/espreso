//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#include "IterSolver.h"

IterSolver::IterSolver():
	timing			("Main CG loop timing "),
	preproc_timing	("Preprocessing timing "),
	postproc_timing	("Postprocessing timing"),
	timeEvalAppa	("Apply Kplus timing "),
	apa_B1t			("x = B1t * lambda "),
	apa_kplus		("multKplus(local or global) "),
	apa_B1			("lambda = B1 * x "),
	apa_allred		("All_Reduce_lambdas "),
	timeEvalPrec	("Apply Precond. timing "),
	prec_kplus		("B1 * P * B1t "),
	prec_allred		("All_Reduce_lambdas "),
	timeEvalProj	("Projector timing "),
	proj_G1t		("x = G1 * lambda "),
	proj_Gthr		("MPI_gather - collective "),
	proj_GGt		("GGt Solve on master node "),
	proj_Sctr		("MPI_Scatter - collective "),
	proj_Gx			("lambda = G1t * x "),
	proj_allred		("All_Reduce_lambdas "),
	ddot_time		("Parallel DDOT - alpha and gamma"),
	proj_time		("Projector_l "),
	appA_time		("ApplyA_l "),
	vec_time		("vector processing in CG "),
	norm_time		("parallel DDOT - norm "),
	proj1_time		("Projector_l - before PREC "),
	proj2_time		("Projector_l - after PREC "),
	prec_time		("Preconditioner "),
	ddot_alpha		("2x ddot for Alpha "),
	ddot_beta		("2x ddot for Beta ")
{

	// Timing objects
	// Main timing object for main CG loop
	timeEvalAppa.addEvent(apa_B1t);
	timeEvalAppa.addEvent(apa_kplus);
	timeEvalAppa.addEvent(apa_B1);
	timeEvalAppa.addEvent(apa_allred);

	timeEvalPrec.addEvent(prec_kplus);
	timeEvalPrec.addEvent(prec_allred);

	timeEvalProj.addEvent(proj_G1t);
	timeEvalProj.addEvent(proj_Gthr);
	timeEvalProj.addEvent(proj_GGt);
	timeEvalProj.addEvent(proj_Sctr);
	timeEvalProj.addEvent(proj_Gx);
	timeEvalProj.addEvent(proj_allred);
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

}

void IterSolver::Preprocessing ( Cluster & cluster )
{

	preproc_timing.totalTime.start();

	// ****************************************************************************
	// *** Coarse problem - Make GGt **********************************************

	TimeEvent createGGT_time("Time to create GGt");
	createGGT_time.start();

	if (USE_DYNAMIC == 0) {
		if (USE_GGtINV == 1)
			CreateGGt_inv_dist( cluster );
		else
			CreateGGt    ( cluster );
	}

	createGGT_time.end();
	createGGT_time.printStatMPI();
	preproc_timing.addEvent(createGGT_time);

	// *** END - Make GGt *********************************************************
	// ****************************************************************************

	preproc_timing.totalTime.end();


}

void IterSolver::Solve_singular ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel )
{

	if ( USE_PIPECG == 1 ) {
		Solve_PipeCG_singular_dom( cluster, in_right_hand_side_primal );
	} else {
		Solve_RegCG_singular_dom ( cluster, in_right_hand_side_primal );
	}

	 postproc_timing.totalTime.start();

	 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
	 timeGetSol.start();
	GetSolution_Primal_singular_parallel( cluster, in_right_hand_side_primal, out_primal_solution_parallel );
	 timeGetSol.endWithBarrier();
	 postproc_timing.addEvent(timeGetSol);

	 postproc_timing.totalTime.endWithBarrier();




}

void IterSolver::Solve_non_singular ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel )
{

	if ( USE_PIPECG == 1 ) {
		Solve_PipeCG_nonsingular ( cluster, in_right_hand_side_primal, out_primal_solution_parallel);
	} else {
		Solve_RegCG_nonsingular  ( cluster, in_right_hand_side_primal, out_primal_solution_parallel);
	}

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

void IterSolver::GetSolution_Primal_singular_parallel  ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
		SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out ) {

	MakeSolution_Primal_singular_parallel(cluster, in_right_hand_side_primal, primal_solution_out );
	//primal_solution_out = primal_solution_parallel;

}

void IterSolver::MakeSolution_Primal_singular_parallel ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
		SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out )  {

	//primal_solution_parallel.clear();
	primal_solution_out.clear();

	// R * mu
	SEQ_VECTOR<SEQ_VECTOR<double> > R_mu_prim_cluster;

	for (eslocal d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR <double > tmp (cluster.domains[d].domain_prim_size);
		if (USE_HFETI == 1)

			if ( esconfig::solver::REGULARIZATION == 0 )
				cluster.domains[d].Kplus_R.DenseMatVec(amplitudes, tmp, 'N', 0, 0);
		  	else
		  		cluster.domains[d].Kplus_Rb.DenseMatVec(amplitudes, tmp, 'N', 0, 0);

		else
			cluster.domains[d].Kplus_R.DenseMatVec(amplitudes, tmp, 'N', d * cluster.domains[d].Kplus_R.cols, 0);

		R_mu_prim_cluster.push_back(tmp);
	}

	for (eslocal d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
		SEQ_VECTOR < double > tmp      ( cluster.domains[d].domain_prim_size  );

		for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = dual_soultion_compressed_parallel[ cluster.domains[d].lambda_map_sub_local[i]];

		cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, tmp, 'T');

		for (eslocal i = 0; i < tmp.size(); i++)
			tmp[i] = in_right_hand_side_primal[d][i] - tmp[i];
			//tmp[i] = cluster.domains[d].f[i] - tmp[i];

		//primal_solution_parallel.push_back(tmp);
		primal_solution_out.push_back(tmp);
	}

	if ( cluster.USE_HFETI == 0) {
		for (eslocal d = 0; d < cluster.domains.size(); d++)
			cluster.domains[d].multKplusLocal(primal_solution_out[d]);
			//cluster.domains[d].multKplusLocal(primal_solution_parallel[d]);
	} else  {
		//cluster.multKplusGlobal_l(primal_solution_parallel);
		cluster.multKplusGlobal_l(primal_solution_out);
	}

 	for (eslocal d = 0; d < cluster.domains.size(); d++) {
		//for (eslocal i = 0; i < primal_solution_parallel[d].size()	; i++) {
 		for (eslocal i = 0; i < primal_solution_out[d].size()	; i++) {
 			primal_solution_out[d][i] = primal_solution_out[d][i] + R_mu_prim_cluster[d][i];
			//primal_solution_parallel[d][i] = primal_solution_parallel[d][i] + R_mu_prim_cluster[d][i];
			//primal_solution_parallel[d][i] = cluster.domains[d].up0[i] + R_mu_prim_cluster[d][i];
		}
	}

}


// *** Singular CG Solvers ***********************************************
void IterSolver::Solve_RegCG_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{

	eslocal dl_size = cluster.my_lamdas_indices.size();

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

	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );

	// *** CG start ***************************************************************

	// t1 = Uc\(Lc\d);
	// x = Ct * t1;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	//double x_norm_l = parallel_norm_compressed(cluster, cluster.vec_d);
	//printf (       "Test probe 1: norm = %1.30f \n", x_norm_l );
	//x_norm_l = parallel_norm_compressed(cluster, x_l);
	//printf (       "Test probe 1: norm = %1.30f \n", x_norm_l );

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	// *** up0 pro ukoncovani v primaru
	////cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
	////	cluster.domains[d].BtLambda_i = cluster.x_prim_cluster2[d];
	////	//cluster.domains[d].BtLambda_i.resize(cluster.domains[d].up0.size(), 0);
	////	cluster.domains[d].norm_f = 0.0;
	////	for (eslocal i = 0; i < cluster.domains[d].up0.size(); i++ ) {
	////		cluster.domains[d].up0[i]     = cluster.domains[d].up0[i] - cluster.x_prim_cluster1[d][i];  // (K+ * f) - (K+ * Bt * lambda)
	////
	////		cluster.domains[d].norm_f += cluster.domains[d].f[i] * cluster.domains[d].f[i];
	////	}
	////}

	// Get norm of f (right hand side)
	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (eslocal d = 0; d < cluster.domains.size(); d++)
		norm_prim_fl += cluster.domains[d].norm_f;

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Reduce   (&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	norm_prim_fg = sqrt(norm_prim_fg);



	// *** r = b - Ax *************************************************************
	cilk_for (eslocal i = 0; i < r_l.size(); i++)
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

		timing.totalTime.start();

		cilk_for (eslocal i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		if (USE_PREC > 0) {

			proj1_time.start();
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj1_time.end();

			// Scale
			prec_time.start();
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, z_l);
			prec_time.end();
			// Re-Scale

			proj2_time.start();
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, z_l, y_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, z_l, y_l, 0 );
			}
			proj2_time.end();

		} else {

			proj_time.start();
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj_time.end();

			cilk_for (eslocal i = 0; i < w_l.size(); i++)
				y_l[i] = w_l[i];

		}


		//------------------------------------------
		if (iter == 0) {									// if outputs.n_it==1;

			cilk_for (eslocal i = 0; i < y_l.size(); i++)
				p_l[i] = y_l[i];							// p = y;

		} else {

			ddot_beta.start();
			beta_l =          parallel_ddot_compressed(cluster, y_l, w_l);
			beta_l = beta_l / parallel_ddot_compressed(cluster, yp_l, wp_l);
			ddot_beta.end();

			cilk_for (eslocal i = 0; i < p_l.size(); i++)
				p_l[i] = y_l[i] + beta_l * p_l[i];			// p = y + beta * p;

		}



		//------------------------------------------
		 appA_time.start();
		apply_A_l_comp_dom_B(timeEvalAppa, cluster, p_l, Ap_l); // apply_A_l_compB(timeEvalAppa, cluster, p_l, Ap_l);
		 appA_time.end();

		//------------------------------------------
		 ddot_alpha.start();
		alpha_l =           parallel_ddot_compressed(cluster, y_l, w_l);
		alpha_l = alpha_l / parallel_ddot_compressed(cluster, p_l, Ap_l);
		 ddot_alpha.end();

		//-----------------------------------------
		// *** up0 pro ukoncovani v primaru

		//// //cilk_
		////for (eslocal d = 0; d < cluster.domains.size(); d++) {
		////	for (eslocal i = 0; i < cluster.domains[d].up0.size(); i++ ) {
		////		cluster.domains[d].up0[i]        -= alpha_l * cluster.x_prim_cluster1[d][i];
		////		cluster.domains[d].BtLambda_i[i] += alpha_l * cluster.x_prim_cluster2[d][i];
		////	}
		////	cluster.domains[d].norm_vec.resize(cluster.domains[d].up0.size());
		////	cluster.domains[d].K.MatVec(cluster.domains[d].up0, cluster.domains[d].norm_vec, 'N');
		////	cluster.domains[d].norm_c = 0.0;
		////	for (eslocal i = 0; i < cluster.domains[d].up0.size(); i++ ) {
		////		cluster.domains[d].norm_vec[i] = cluster.domains[d].norm_vec[i]
		////			                           + cluster.domains[d].BtLambda_i[i]
		////									   - cluster.domains[d].f[i];
		////
		////		cluster.domains[d].norm_c += cluster.domains[d].norm_vec[i] * cluster.domains[d].norm_vec[i];
		////	}
		////}

		//double norm_prim_l = 0.0;
		//double norm_prim_g = 0.0;
		//for (eslocal d = 0; d < cluster.domains.size(); d++)
		//	norm_prim_l += cluster.domains[d].norm_c;

		//MPI_Allreduce(&norm_prim_l, &norm_prim_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		////MPI_Reduce(&norm_prim_l, &norm_prim_g, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		//norm_prim_g = sqrt(norm_prim_g);





		//------------------------------------------
		cilk_for (eslocal i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * Ap_l[i];
		}

		timing.totalTime.end();
		timing.totalTime.printLastStatMPI();

		 norm_time.start();
		norm_l = parallel_norm_compressed(cluster, w_l);
		 norm_time.end();


		if (mpi_rank == mpi_root) {
			//printf (       "Iter MPI %d - norm dual %1.20f - tol %1.20f - norm prim %1.20f norm f %1.20f : %1.20f \n", iter+1, norm_l, tol, norm_prim_g, norm_prim_fg, norm_prim_g / norm_prim_fg);
//		printf (       "Iter MPI %5d - norm dual %1.20f - tol %1.20f \n", iter+1, norm_l, tol );
//



      int my_prec= 10; //log10(int(1./epsilon));
      std::cout.clear();
      std::cout<<"Iter MPI ";
      std::cout<<std::setw(my_prec+2);
      std::cout<<iter+1;
      std::cout.precision(my_prec+2);
      std::cout<<" ||normed_residual|| = "<<norm_l/tol*epsilon;
      std::cout<<" ||residual|| = "<<norm_l;
      std::cout.precision(my_prec);
      std::cout<<",   epsilon = "<<epsilon <<"\n";
      //printf (       "Iter MPI %5d - ||normed_residual|| %1.20f - tol %1.20f\n", iter+1, norm_l/tol*epsilon,epsilon );
    //  std::cout<<      "Iter MPI "<< iter+1 << "- norm dual "<< norm_l << " - tol " <<tol <<"\n";

			//if (log_active == 1)
			//fprintf(stream,"Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
		}
//    std::cout.setstate(std::ios_base::failbit);


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


	// *** Preslocal out the timing for the iteration loop ***************************************

	if (USE_PREC > 0) {
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
	} else {
		timing.addEvent(proj_time);
	}

	timing.addEvent(appA_time );
	timing.addEvent(ddot_beta);
	timing.addEvent(ddot_alpha);

	// *** END - Preslocal out the timing for the iteration loop ***********************************

}

void IterSolver::Solve_PipeCG_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
	eslocal dl_size = cluster.my_lamdas_indices.size();

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

	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );

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
	if (USE_PREC > 0) {

		cilk_for (eslocal i = 0; i < r_l.size(); i++)
			tmp_l[i] = b_l[i] - Ax_l[i];

		if (USE_GGtINV == 1)
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, r_l, 0 );
		else
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, r_l, 0 );

	} else {

		cilk_for (eslocal i = 0; i < r_l.size(); i++)
			r_l[i] = b_l[i] - Ax_l[i];
	}


	if (USE_PREC > 0) {

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


	if (USE_PREC > 0) {
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

		timing.totalTime.start();

		alpha_lp = alpha_l;
		gama_lp  = gama_l;

		//------------------------------------------
		ddot_time.start();
		MPI_Request mpi_req;
		MPI_Status mpi_stat;

		SEQ_VECTOR <double> reduction_tmp (2,0);
		SEQ_VECTOR <double> send_buf (2,0);
		parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, &mpi_req, reduction_tmp, send_buf);

		ddot_time.end();

		if (USE_PREC > 0) {

			prec_time.start();
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, tmp_l);
			prec_time.end();

			proj_time.start();
			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, m_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, m_l, 0 );
			proj_time.end();

			appA_time.start();
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, tmp_l); //apply_A_l_compB(timeEvalAppa, cluster, m_l, tmp_l);
			appA_time.end();

			proj_time.start();
			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, n_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, n_l, 0 );
			proj_time.end();

		} else {

			//------------------------------------------
			proj_time.start();

			if (USE_GGtINV == 1)
				Projector_l_inv_compG( timeEvalProj, cluster, w_l, m_l, 0 );
			else
				Projector_l_compG    ( timeEvalProj, cluster, w_l, m_l, 0 );

			proj_time.end();

			//------------------------------------------
			appA_time.start();
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, n_l); //apply_A_l_compB(timeEvalAppa, cluster, m_l, n_l);
			appA_time.end();
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
		vec_time.start();
		if (iter == 0) {
			beta_l  = 0;
			alpha_l = gama_l / delta_l;
		} else {
			beta_l  = gama_l / gama_lp;
			alpha_l = gama_l / (delta_l - beta_l * gama_l / alpha_lp);
		}

		cilk_for (eslocal i = 0; i < r_l.size(); i++) {
			z_l[i] = n_l[i] + beta_l  * z_l[i];
			q_l[i] = m_l[i] + beta_l  * q_l[i];
			s_l[i] = w_l[i] + beta_l  * s_l[i];
			p_l[i] = u_l[i] + beta_l  * p_l[i];
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * s_l[i];
			u_l[i] = u_l[i] - alpha_l * q_l[i];
			w_l[i] = w_l[i] - alpha_l * z_l[i];
		}
		vec_time.end();


		timing.totalTime.end();
		timing.totalTime.printLastStatMPI();

		norm_time.start();

		// POZOR - tady se to ukoncuje jinak = musime probrat
		if (USE_PREC > 0)
			norm_l = parallel_norm_compressed(cluster, r_l);
		else
			norm_l = parallel_norm_compressed(cluster, u_l);

		norm_time.end();

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

	cilk_for(eslocal i = 0; i < r_l.size(); i++)
		r_l[i] = b_l[i] - Ax_l[i];

	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************
	timing.addEvent(ddot_time);

	// *** Preslocal out the timing for the iteration loop ***************************************

	if (USE_PREC > 0) {
		//timing.addEvent(proj1_time);
		timing.addEvent(proj_time);
		timing.addEvent(prec_time );
		//timing.addEvent(proj2_time);
	} else {
		timing.addEvent(proj_time);
	}

	timing.addEvent(appA_time);
	timing.addEvent(vec_time );

}


// *** Non-singular CG Solvers *******************************************
void IterSolver::Solve_RegCG_nonsingular  ( Cluster & cluster,
										    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
										    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel) {

	eslocal dl_size = cluster.my_lamdas_indices.size();

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
//	for (eslocal d = 1; d < cluster.domains.size(); d++) {
//		cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[0]);
//		cluster.domains[d].B1_comp.MatVec(cluster.x_prim_cluster1[0], cluster.compressed_tmp, 'N', 0, 0, 1.0);
//	}
//	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);

	//// *** convert right hand side to dual
	std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
	SEQ_VECTOR < double > y_out_tmp;
	for (eslocal d = 0; d < cluster.domains.size(); d++) {
		// *** convert right hand side to dual
		cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[d]);

		y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
		cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

		for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
	}
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);



	int iter = 0;
	timing.totalTime.reset();

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);

	cilk_for (eslocal i = 0; i < r_l.size(); i++) {	// r = b - Ax;
		r_l[i] = b_l[i] - Ax_l[i];
		wp_l[i] = 0.0;
		yp_l[i] = 0.0;
	}

	tol = epsilon * parallel_norm_compressed(cluster, b_l);

	for ( iter = 0; iter < 1000; iter++) {
		timing.totalTime.start();

		cilk_for (eslocal i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		if (USE_PREC > 0) {
			cilk_for (eslocal i = 0; i < w_l.size(); i++)
				w_l[i] = r_l[i];
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, y_l);
		} else {
			cilk_for (eslocal i = 0; i < w_l.size(); i++)
				w_l[i] = r_l[i];
			cilk_for (eslocal i = 0; i < w_l.size(); i++)
				y_l[i] = w_l[i];
		}


		if (iter == 0) {									// if outputs.n_it==1;
			cilk_for (eslocal i = 0; i < y_l.size(); i++)
				p_l[i] = y_l[i];							// p = y;
		} else {
			ddot_beta.start();
			beta_l =          parallel_ddot_compressed(cluster, y_l, w_l);
			beta_l = beta_l / parallel_ddot_compressed(cluster, yp_l, wp_l);
			ddot_beta.end();

			cilk_for (eslocal i = 0; i < p_l.size(); i++)
				p_l[i] = y_l[i] + beta_l * p_l[i];			// p = y + beta * p;
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, p_l, Ap_l);



		alpha_l =           parallel_ddot_compressed(cluster, y_l, w_l);
		alpha_l = alpha_l / parallel_ddot_compressed(cluster, p_l, Ap_l);

		cilk_for (eslocal i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * Ap_l[i];
		}

		timing.totalTime.end();

		norm_l = parallel_norm_compressed(cluster, r_l);

		if (mpi_rank == mpi_root) {
			printf ( "CG Iter %d - norm %1.30f - tol %1.30f \r" , iter+1, norm_l, tol);
			//if (log_active == 1)
			//	fprintf(stream,"Iter MPI %d - norm %1.30f - tol %1.30f \n", iter+1, norm_l, tol);
		}


		if (norm_l < tol)
			break;

	} // end iter loop

	if (mpi_rank == mpi_root) cout << endl;

	 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
	 timeGetSol.start();

	// reconstruction of u
	cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );

		for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_l[ cluster.domains[d].lambda_map_sub_local[i]];

		cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');

		for(eslocal i = 0; i < in_right_hand_side_primal[d].size(); i++)
			out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];

		cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);

	}

	 timeGetSol.endWithBarrier();
	 timing.addEvent(timeGetSol);


}

void IterSolver::Solve_PipeCG_nonsingular ( Cluster & cluster,
											SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
											SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel) {

		eslocal dl_size = cluster.my_lamdas_indices.size();

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
//		for (eslocal d = 1; d < cluster.domains.size(); d++) {
//			cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[0]);
//			cluster.domains[d].B1_comp.MatVec(cluster.x_prim_cluster1[0], cluster.compressed_tmp, 'N', 0, 0, 1.0);
//		}
//		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);

		//// *** convert right hand side to dual
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			// *** convert right hand side to dual
			cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[d]);

			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);


		int iter = 0;
		timing.totalTime.reset();

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);

		cilk_for (eslocal i = 0; i < r_l.size(); i++) {	// r = b - Ax;
			r_l[i] = b_l[i] - Ax_l[i];
			wp_l[i] = 0.0;
			yp_l[i] = 0.0;
		}

		if (USE_PREC > 0) {
			apply_prec_comp_dom_B(timeEvalPrec, cluster, r_l, u_l);
		} else {
			cilk_for (eslocal i = 0; i < r_l.size(); i++)
				u_l = r_l;
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, w_l);

		tol = epsilon * parallel_norm_compressed(cluster, b_l);

		for ( iter = 0; iter < CG_max_iter; iter++) {

			timing.totalTime.start();

			alpha_lp = alpha_l;
			gama_lp  = gama_l;

			ddot_time.start();
			MPI_Request mpi_req;
			MPI_Status mpi_stat;

			SEQ_VECTOR <double> reduction_tmp (2,0);
			SEQ_VECTOR <double> send_buf (2,0);
			parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, &mpi_req, reduction_tmp, send_buf);

			ddot_time.end();

			if (USE_PREC > 0) {

				prec_time.start();
				apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, m_l);
				prec_time.end();

			} else {

				cilk_for (eslocal i = 0; i < m_l.size(); i++)
					m_l[i] = w_l[i];

			}

			//------------------------------------------
			appA_time.start();
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, n_l);
			appA_time.end();
			//------------------------------------------


#ifndef WIN32
#ifdef USE_MPI_3
			MPI_Wait(&mpi_req, &mpi_stat);
#endif
#endif
			gama_l  = reduction_tmp[0];
			delta_l = reduction_tmp[1];

			//------------------------------------------
			vec_time.start();
			if (iter == 0) {
				beta_l  = 0;
				alpha_l = gama_l / delta_l;
			} else {
				beta_l  = gama_l / gama_lp;
				alpha_l = gama_l / (delta_l - beta_l * gama_l / alpha_lp);
			}

			cilk_for (eslocal i = 0; i < r_l.size(); i++) {
				z_l[i] = n_l[i] + beta_l  * z_l[i];
				q_l[i] = m_l[i] + beta_l  * q_l[i];
				s_l[i] = w_l[i] + beta_l  * s_l[i];
				p_l[i] = u_l[i] + beta_l  * p_l[i];
				x_l[i] = x_l[i] + alpha_l * p_l[i];
				r_l[i] = r_l[i] - alpha_l * s_l[i];
				u_l[i] = u_l[i] - alpha_l * q_l[i];
				w_l[i] = w_l[i] - alpha_l * z_l[i];
			}
			vec_time.end();


			timing.totalTime.end();
			//timing.totalTime.printLastStatMPI();

			norm_time.start();
			//norm_l = parallel_norm_compressed(cluster, u_l);
			norm_l = parallel_norm_compressed(cluster, r_l);
			norm_time.end();

			if (mpi_rank == mpi_root) {
				//printf (       "Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
				printf ( "CG Iter %d - norm %1.30f - tol %1.30f \r" , iter+1, norm_l, tol);
				//if (log_active == 1)
				//fprintf(stream,"Iter MPI %d - norm %1.20f - tol %1.20f \n", iter+1, norm_l, tol);
			}

			if (norm_l < tol)
				break;

		} // end iter loop

		// reconstruction of u

		if (mpi_rank == mpi_root) cout << endl;

		 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
		 timeGetSol.start();

		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );

			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_l[ cluster.domains[d].lambda_map_sub_local[i]];

			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');

			for(eslocal i = 0; i < in_right_hand_side_primal[d].size(); i++)
				out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];

			cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);

		}

		 timeGetSol.endWithBarrier();
		 timing.addEvent(timeGetSol);


}


// *** Coarse problem routines *******************************************
void IterSolver::CreateGGt( Cluster & cluster )

{

	double start = omp_get_wtime();

	double sc1 = omp_get_wtime();
	SparseMatrix G;


	//if (mpi_rank == mpi_root)
	//	G.MatAppend(cluster.G1);

	//for (eslocal mr = 1; mr < mpi_size; mr++) {
	//	SparseMatrix Gtmp;
	//	SendMatrix(mpi_rank, mr, cluster.G1, mpi_root, Gtmp);

	//	if (mpi_rank == mpi_root) {
	//		G.MatAppend(Gtmp);
	//		Gtmp.Clear();
	//	}
	//}

	//// **** Log N MPI reduce
	int count_cv = 0;
	for (eslocal li = 2; li <= 2*mpi_size; li = li * 2 ) {

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
			//printf(" Collecting matrices G : %d of %d \r", count_cv, mpi_size);
			std::cout<<" Collecting matrices G : " << count_cv <<" of " << mpi_size << "\n";

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

		MKL_Set_Num_Threads(PAR_NUM_THREADS);
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
		std::stringstream ss;
		ss << "Create GGt -> rank: " << esconfig::MPIrank;
		GGt.Factorization(ss.str());
		cout << "Factorization = " << omp_get_wtime() - t1 << endl;


		t1 = omp_get_wtime();
		GGt.msglvl = 0;
		//TODO:
		MKL_Set_Num_Threads(1);
	}

	double ep1 = omp_get_wtime();


	if (mpi_rank == mpi_root)
		GGtsize = GGt.cols;


	MPI_Bcast( & GGtsize, 1, esglobal_mpi, 0, MPI_COMM_WORLD);


#if TIME_MEAS >= 1
	eslocal rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double end = omp_get_wtime();
		cout <<	"CG Loop - Create GGt  - collect all matrices   - Runtime = " << ec1 - sc1 << " s \n";
		cout <<	"CG Loop - Create GGt  - GGt fact. processing   - Runtime = " << ep1 - sp1 << " s \n";
		cout <<	"CG Loop - Create GGt  - total = proc + comm    - Runtime = " << end - start << " s \n\n";
	}
#endif

}

void IterSolver::CreateGGt_inv_dist( Cluster & cluster )
{

	// temp variables
	vector < SparseMatrix > G_neighs   ( cluster.my_neighs.size() );
	vector < SparseMatrix > GGt_neighs ( cluster.my_neighs.size() );
	SparseMatrix Gt_l;
	SparseMatrix GGt_l;
	SparseMatrix GGt_Mat_tmp;
	SparseSolverCPU GGt_tmp;


        /* Numbers of processors, value of OMP_NUM_THREADS */
        int num_procs;
        char * var = getenv("PAR_NUM_THREADS");
        if(var != NULL)
                sscanf( var, "%d", &num_procs );
        else {
                //printf("Set environment PAR_NUM_THREADS to 1");
          std::cerr<<"Set environment PAR_NUM_THREADS to 1\n";
                exit(1);
        }
        GGt_tmp.iparm[2]  = num_procs;


	 TimeEvent SaRGlocal("Send a Receive local G1 matrices to neighs. "); SaRGlocal.start();
	for (eslocal neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ )
		SendMatrix(cluster.G1, cluster.my_neighs[neigh_i]);

	for (eslocal neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ )
		RecvMatrix(G_neighs[neigh_i], cluster.my_neighs[neigh_i]);
	 SaRGlocal.end(); SaRGlocal.printStatMPI(); preproc_timing.addEvent(SaRGlocal);


	 TimeEvent Gt_l_trans("Local G1 matrix transpose to create Gt "); Gt_l_trans.start();
	if (cluster.USE_HFETI == 0)
		cluster.G1.MatTranspose(Gt_l);
	 Gt_l_trans.end(); Gt_l_trans.printStatMPI(); preproc_timing.addEvent(Gt_l_trans);

	 TimeEvent GxGtMatMat("Local G x Gt MatMat "); GxGtMatMat.start();

	if (cluster.USE_HFETI == 0)
		GGt_l.MatMat(cluster.G1, 'N', Gt_l);
	else
		GGt_l.MatMatT(cluster.G1, cluster.G1);

	 GxGtMatMat.end(); GxGtMatMat.printStatMPI(); preproc_timing.addEvent(GxGtMatMat);
	 //GxGtMatMat.PrintLastStatMPI_PerNode(0.0);


	for (eslocal i = 0; i < GGt_l.CSR_J_col_indices.size(); i++)
		GGt_l.CSR_J_col_indices[i] += mpi_rank * cluster.G1.rows;

	GGt_l.cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;

	 TimeEvent GGTNeighTime("G1t_local x G1_neigh MatMat(N-times) "); GGTNeighTime.start();
	cilk_for (eslocal neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {

		if (cluster.USE_HFETI == 0)
			GGt_neighs[neigh_i].MatMat(G_neighs[neigh_i], 'N', Gt_l);
		else
			GGt_neighs[neigh_i].MatMatT(G_neighs[neigh_i], cluster.G1);

		GGt_neighs[neigh_i].MatTranspose();

		eslocal inc = cluster.G1.rows * cluster.my_neighs[neigh_i];
		for (eslocal i = 0; i < GGt_neighs[neigh_i].CSR_J_col_indices.size(); i++)
			GGt_neighs[neigh_i].CSR_J_col_indices[i] += inc;

		GGt_neighs[neigh_i].cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;
		G_neighs[neigh_i].Clear();
	}
	 GGTNeighTime.end(); GGTNeighTime.printStatMPI(); preproc_timing.addEvent(GGTNeighTime);
	 //GGTNeighTime.PrintLastStatMPI_PerNode(0.0);

	 TimeEvent GGtLocAsm("Assembling row of GGt per node - MatAddInPlace "); GGtLocAsm.start();
	for (eslocal neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		GGt_l.MatAddInPlace(GGt_neighs[neigh_i], 'N', 1.0);
		GGt_neighs[neigh_i].Clear();
	}
	 GGtLocAsm.end(); GGtLocAsm.printStatMPI(); preproc_timing.addEvent(GGtLocAsm);


	TimeEvent collectGGt_time("Collect GGt pieces to master"); 	collectGGt_time.start();
	int count_cv_l = 0;
	for (eslocal li = 2; li <= 2*mpi_size; li = li * 2 ) {

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
//			printf(" Collecting matrices G : %d of %d \r", count_cv_l, mpi_size);
			std::cout<<" Collecting matrices G : " << count_cv_l <<" of " << mpi_size << "\n";

	}
	collectGGt_time.end(); collectGGt_time.printStatMPI(); preproc_timing.addEvent(collectGGt_time);

	if (mpi_rank == 0) cout << endl;

	if (mpi_rank == 0)  {
		GGt_Mat_tmp.RemoveLower();
		SpyText(GGt_Mat_tmp);
	}

	MKL_Set_Num_Threads(PAR_NUM_THREADS);

	TimeEvent GGt_bcast_time("Time to broadcast GGt from master all"); GGt_bcast_time.start();
	BcastMatrix(mpi_rank, mpi_root, mpi_root, GGt_Mat_tmp);
	GGt_bcast_time.end(); GGt_bcast_time.printStatMPI(); preproc_timing.addEvent(GGt_bcast_time);

	// Create Sparse Direct solver for GGt
	if (mpi_rank == mpi_root)
		GGt_tmp.msglvl = 1;

	TimeEvent importGGt_time("Time to import GGt matrix into solver"); importGGt_time.start();
	GGt_tmp.ImportMatrix(GGt_Mat_tmp);
	importGGt_time.end(); importGGt_time.printStatMPI(); preproc_timing.addEvent(importGGt_time);

	GGt_Mat_tmp.Clear();

	TimeEvent GGtFactor_time("GGT Factorization time"); GGtFactor_time.start();
	GGt_tmp.SetThreaded();
	std::stringstream ss;
	ss << "Create GGt_inv_dist-> rank: " << esconfig::MPIrank;
	GGt_tmp.Factorization(ss.str());
	GGtFactor_time.end(); GGtFactor_time.printStatMPI(); preproc_timing.addEvent(GGtFactor_time);

	TimeEvent GGT_rhs_time("Time to create RHS for get GGTINV"); GGT_rhs_time.start();
	SEQ_VECTOR <double> rhs   (cluster.G1.rows * GGt_tmp.rows, 0);
	cluster.GGtinvV.resize(cluster.G1.rows * GGt_tmp.rows, 0);

	for (eslocal i = 0; i < cluster.G1.rows; i++) {
		eslocal index = (GGt_tmp.rows * i) + (cluster.G1.rows * mpi_rank) + i;
		rhs[index] = 1;
	}
	GGT_rhs_time.end(); GGT_rhs_time.printStatMPI(); preproc_timing.addEvent(GGT_rhs_time);

	TimeEvent GGt_solve_time("Running solve to get stripe(s) of GGtINV"); GGt_solve_time.start();

	GGt_tmp.Solve(rhs, cluster.GGtinvV, cluster.G1.rows);

	cluster.GGtinvM.dense_values = cluster.GGtinvV;
	cluster.GGtinvM.cols = cluster.G1.rows;
	cluster.GGtinvM.rows = GGt_tmp.rows;
        cluster.GGtinvM.type = 'G';

	GGtsize  = GGt_tmp.cols;

	GGt.cols = GGt_tmp.cols;
	GGt.rows = GGt_tmp.rows;
	GGt.nnz  = GGt_tmp.nnz;

	GGt_tmp.msglvl = 0;
	GGt_tmp.Clear();

	GGt_solve_time.end(); GGt_solve_time.printStatMPI(); preproc_timing.addEvent(GGt_solve_time);

	MKL_Set_Num_Threads(1);

}


// *** Projector routines ************************************************
void IterSolver::Projector_l_compG (TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
{

	time_eval.totalTime.start();

	//eslocal dual_size    = cluster.domains[0].B1.rows;
	eslocal d_local_size = cluster.G1_comp.rows;
	eslocal mpi_root     = 0;

	SEQ_VECTOR<double> d_local( d_local_size );
	SEQ_VECTOR<double> d_mpi  ( GGtsize );

	time_eval.timeEvents[0].start();
	if ( output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 1)
		d_local = x_in;
	else
		cluster.G1_comp.MatVec(x_in, d_local, 'N');
	time_eval.timeEvents[0].end();



	time_eval.timeEvents[1].start();
	MPI_Gather(&d_local[0], d_local_size, MPI_DOUBLE,
		&d_mpi[0], d_local_size, MPI_DOUBLE,
		mpi_root, MPI_COMM_WORLD);
	time_eval.timeEvents[1].end();


//	for (int i = 0; i < d_mpi.size(); i++)
//	printf (       "Test probe 1: %d norm = %1.30f \n", i, d_mpi[i] );


	time_eval.timeEvents[2].start();
	if (mpi_rank == mpi_root ) {
		GGt.Solve(d_mpi);				// t1 = Uc\(Lc\d);
	}
	time_eval.timeEvents[2].end();

//	for (int i = 0; i < d_mpi.size(); i++)
//	printf (       "Test probe 1: %d norm = %1.30f \n", i, d_mpi[i] );


	time_eval.timeEvents[3].start();
	MPI_Scatter( &d_mpi[0],      d_local_size, MPI_DOUBLE,
		&d_local[0], d_local_size, MPI_DOUBLE,
		mpi_root, MPI_COMM_WORLD);
	time_eval.timeEvents[3].end();

	if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 2) {
		// for mu calculation
		y_out = d_local;

	} else {

		time_eval.timeEvents[4].start();
		//cluster.G1t_comp.MatVec(d_local, cluster.compressed_tmp, 'N'); // SUPER POZOR
		cluster.G1_comp.MatVec(d_local, cluster.compressed_tmp, 'T');
		time_eval.timeEvents[4].end();

		time_eval.timeEvents[5].start();
		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
		time_eval.timeEvents[5].end();

		if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 0) {
			cilk_for (eslocal i = 0; i < x_in.size(); i++)
				y_out[i] = x_in[i] - y_out[i];
		}

	}

	time_eval.totalTime.end();

}

void IterSolver::Projector_l_inv_compG (TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
{

	time_eval.totalTime.start();

	eslocal d_local_size = cluster.G1_comp.rows;
	eslocal mpi_root     = 0;

	SEQ_VECTOR<double> d_local( d_local_size );
	SEQ_VECTOR<double> d_mpi  ( GGtsize );

	time_eval.timeEvents[0].start();
	if ( output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 1)
		d_local = x_in;
	else
		cluster.G1_comp.MatVec(x_in, d_local, 'N');
	time_eval.timeEvents[0].end();

	time_eval.timeEvents[1].start();
	MPI_Allgather(&d_local[0], d_local_size, MPI_DOUBLE,
		&d_mpi[0], d_local_size, MPI_DOUBLE,
		MPI_COMM_WORLD);
	time_eval.timeEvents[1].end();

	time_eval.timeEvents[2].start();
	//for (eslocal j = 0; j < cluster.G1_comp.rows; j++) {
	//	d_local[j] = 0;
	//	for (eslocal i = 0; i < GGt.rows; i++ ) {
	//		d_local[j] += cluster.GGtinvV[j * GGt.rows + i] * d_mpi[i];			// t1 = Uc\(Lc\d);
	//	}
	//}
	cluster.GGtinvM.DenseMatVec(d_mpi, d_local, 'T');
	time_eval.timeEvents[2].end();

	time_eval.timeEvents[3].start();
	//MPI_Scatter( &d_mpi[0],      d_local_size, MPI_DOUBLE,
	//	&d_local[0], d_local_size, MPI_DOUBLE,
	//	mpi_root, MPI_COMM_WORLD);
	time_eval.timeEvents[3].end();

	if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 2) {
		y_out = d_local; // for RBM amplitudes calculation
	} else {

		time_eval.timeEvents[4].start();
		//cluster.G1t_comp.MatVec(d_local, cluster.compressed_tmp, 'N'); SUPER POZOR
		cluster.G1_comp.MatVec(d_local, cluster.compressed_tmp, 'T');
		time_eval.timeEvents[4].end();

		time_eval.timeEvents[5].start();
		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
		time_eval.timeEvents[5].end();

		if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 0) {
			cilk_for (eslocal i = 0; i < y_out.size(); i++)
				y_out[i] = x_in[i] - y_out[i];
		}

	}

	time_eval.totalTime.end();

}

// *** Action of K+ routines *********************************************

void IterSolver::apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {
       time_eval.totalTime.start();

       // number of Xeon Phi devices (0 for fallback to CPU)
       eslocal numDevices = cluster.NUM_MICS;

       #ifdef MIC
       if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {

               // HFETI on MIC using Schur

               time_eval.timeEvents[0].start();

               eslocal maxDevNumber = numDevices;
               if (numDevices == 0) {
                       maxDevNumber = 1;               // if number of MICs is zero, we use a CPU
               }
               eslocal matrixPerPack = cluster.domains.size() / maxDevNumber;
               eslocal offset = 0;

               for ( eslocal i = 0; i < maxDevNumber; i++ ) {
                       if ( i == maxDevNumber - 1 ) {
                               // add the remaining domains to the last pack
                               matrixPerPack += cluster.domains.size() % maxDevNumber;
                       }
                       cilk_for (eslocal d = offset; d < offset + matrixPerPack; d++) {
                               SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
                               for ( eslocal j = 0; j < cluster.domains[d].lambda_map_sub_local.size(); j++ ) {
                                       cluster.B1KplusPacks[i].SetX(d - offset, j, x_in[ cluster.domains[d].lambda_map_sub_local[j]]);
                                       x_in_tmp[j] = x_in[ cluster.domains[d].lambda_map_sub_local[j]];
                               }
                               cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
                       }
                       offset += matrixPerPack;
               }

       #pragma omp parallel num_threads( maxDevNumber )
       {
               if ( numDevices > 0 ) {
                       // start async. computation on MICs
                       cluster.B1KplusPacks[ omp_get_thread_num() ].DenseMatsVecsMIC_Start( 'N' );
               } else {
                       // perform matrix-vector multiplication on CPU using MIC data structures
                       cluster.B1KplusPacks[ 0 ].DenseMatsVecs( 'N' );
               }
       }
       time_eval.timeEvents[0].end();

       // perform simultaneous computation on CPU
       time_eval.timeEvents[1].start();
       cluster.multKplusGlobal_Kinv( cluster.x_prim_cluster1 );
       time_eval.timeEvents[1].end();

       time_eval.timeEvents[2].start();
       std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);

       #pragma omp parallel num_threads( maxDevNumber )
       {
               if ( numDevices > 0 ) {
                       // synchronize computation
                       cluster.B1KplusPacks[ omp_get_thread_num() ].DenseMatsVecsMIC_Sync(  );
               }
       }

       // extract the result
       offset = 0;
       matrixPerPack = cluster.domains.size() / maxDevNumber;
       for ( eslocal i = 0; i < maxDevNumber; i++ ) {
               if ( i == maxDevNumber - 1 ) {
                       matrixPerPack += cluster.domains.size() % maxDevNumber;
               }
               cilk_for ( eslocal d = offset ; d < offset + matrixPerPack; d++ ) {
                       cluster.B1KplusPacks[i].GetY(d - offset, cluster.domains[d].compressed_tmp);
               }
               offset+=matrixPerPack;
       }


               time_eval.timeEvents[2].start();
               std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
               for (eslocal d = 0; d < cluster.domains.size(); d++) {
                       for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                               cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
               }
               time_eval.timeEvents[2].end();
               //std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
               SEQ_VECTOR < double > y_out_tmp;
               for (eslocal d = 0; d < cluster.domains.size(); d++) {
                       y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
                       cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

                       for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                               cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
               }
               time_eval.timeEvents[2].end();
       }
       #else


	if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {
		time_eval.timeEvents[0].start();
		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
				x_in_tmp[i]                            = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
#ifdef CUDA
				if (cluster.domains[d].isOnACC == 1) {
					cluster.domains[d].cuda_pinned_buff[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
					//cluster.domains[d].cuda_pinned_buff_fl[i] = (float) x_in[ cluster.domains[d].lambda_map_sub_local[i]];
				}
#endif
			}

#ifdef CUDA
			if (cluster.domains[d].isOnACC == 1) {
				cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start( cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff,'N',0 );
			//cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start_fl( cluster.domains[d].cuda_pinned_buff_fl, cluster.domains[d].cuda_pinned_buff_fl,'N',0 );
			} else {
				cluster.domains[d].B1Kplus.DenseMatVec (x_in_tmp, cluster.domains[d].compressed_tmp);
			}
#else
			cluster.domains[d].B1Kplus.DenseMatVec (x_in_tmp, cluster.domains[d].compressed_tmp);
#endif

			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
		}
		time_eval.timeEvents[0].end();

		time_eval.timeEvents[1].start();
		cluster.multKplusGlobal_Kinv( cluster.x_prim_cluster1 );
		time_eval.timeEvents[1].end();

		time_eval.timeEvents[2].start();
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);

#ifdef CUDA
		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d].isOnACC == 1) {
				cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_sync ( );
			}
		}

		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d].isOnACC == 1) {
				for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
					cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].cuda_pinned_buff[i];
					//cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += (double)cluster.domains[d].cuda_pinned_buff_fl[i];
				}
			} else {
				for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
					cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
				}
			}
		}

#else
		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
		}
#endif

		//std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		time_eval.timeEvents[2].end();

	}
	time_eval.totalTime.start();
	#endif

#ifdef MIC
if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 0) {

       // classical FETI on MIC using Schur

       time_eval.timeEvents[0].start();
       time_eval.timeEvents[0].end();

       time_eval.timeEvents[1].start();

       eslocal maxDevNumber = numDevices;
       if (numDevices == 0) {
               maxDevNumber = 1;       // if number of MICs is zero, we use a CPU
       }

       eslocal matrixPerPack = cluster.domains.size() / maxDevNumber;
       eslocal offset = 0;

       for ( eslocal i = 0; i < maxDevNumber; i++ ) {
               if ( i == maxDevNumber - 1 ) {
                       // add the remaining domains to the last pack
                       matrixPerPack += cluster.domains.size() % maxDevNumber;
               }
               cilk_for (eslocal d = offset; d < offset + matrixPerPack; d++) {
                       //SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
                       for ( eslocal j = 0; j < cluster.domains[d].lambda_map_sub_local.size(); j++ )
                               cluster.B1KplusPacks[i].SetX(d - offset, j, x_in[ cluster.domains[d].lambda_map_sub_local[j]]);
               }
               offset += matrixPerPack;
      }

#pragma omp parallel num_threads( maxDevNumber )
{
       eslocal device = omp_get_thread_num();
       if ( numDevices > 0 ) {
               // run matrix-vector multiplication of MIC
               cluster.B1KplusPacks[ device ].DenseMatsVecsMIC( 'N' );
       } else {
               // run matrix-vector multiplication on CPU
               cluster.B1KplusPacks[ 0 ].DenseMatsVecs( 'N' );
       }
}
       offset = 0;
       matrixPerPack = cluster.domains.size() / maxDevNumber;
       for ( eslocal i = 0; i < maxDevNumber; i++ ) {
               if ( i == maxDevNumber - 1 ) {
                       matrixPerPack += cluster.domains.size() % maxDevNumber;
               }
               cilk_for ( eslocal d = offset ; d < offset + matrixPerPack; d++ ) {
                       cluster.B1KplusPacks[i].GetY(d - offset, cluster.domains[d].compressed_tmp);
               }
               offset+=matrixPerPack;
       }
       time_eval.timeEvents[1].end();

       time_eval.timeEvents[2].start();
       std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
       for (eslocal d = 0; d < cluster.domains.size(); d++) {
               for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                       cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
       }
       time_eval.timeEvents[2].end();
}
#else


	if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 0) {
		time_eval.timeEvents[0].start();
		//// POZOR - jen pro porovnani vykonu CPU a GPU
		//cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
		//	SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
		//	for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
		//		x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
		//	cluster.domains[d].B1Kplus.DenseMatVec(x_in_tmp, cluster.domains[d].compressed_tmp);
		//}
		//// POZOR - jen pro porovnani vykonu CPU a GPU
		time_eval.timeEvents[0].end();

		time_eval.timeEvents[1].start();
		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
#ifdef CUDA
			if (cluster.domains[d].isOnACC == 1) {
				cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy(x_in_tmp, cluster.domains[d].compressed_tmp,'N',0);
			} else {
				cluster.domains[d].B1Kplus.DenseMatVec            (x_in_tmp, cluster.domains[d].compressed_tmp);
			}
#elif MICEXP
			cluster.domains[d].B1Kplus.DenseMatVecMIC_wo_Copy (x_in_tmp, cluster.domains[d].compressed_tmp,'N',0);
#else
			cluster.domains[d].B1Kplus.DenseMatVec            (x_in_tmp, cluster.domains[d].compressed_tmp);
#endif
		}
		time_eval.timeEvents[1].end();

		time_eval.timeEvents[2].start();
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
		}
		time_eval.timeEvents[2].end();
	}

#endif

	if (cluster.USE_KINV == 0) {
		time_eval.timeEvents[0].start();
		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
			//cluster.x_prim_cluster2[d] = cluster.x_prim_cluster1[d]; // POZOR zbytecne kopirovani // prim norm
		}
		time_eval.timeEvents[0].end();

		time_eval.timeEvents[1].start();
		if (cluster.USE_HFETI == 0) {
			cilk_for (eslocal d = 0; d < cluster.domains.size(); d++)
				cluster.domains[d].multKplusLocal(cluster.x_prim_cluster1[d]);
			} else {
				cluster.multKplusGlobal_l(cluster.x_prim_cluster1);
		}
		time_eval.timeEvents[1].end();


		time_eval.timeEvents[2].start();
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		time_eval.timeEvents[2].end();

	}

	time_eval.timeEvents[3].start();
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[3].end();

	time_eval.totalTime.end();

}

void IterSolver::apply_prec_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ) {

	time_eval.totalTime.start();

	time_eval.timeEvents[0].start();

	cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
		for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling

		// Lumped Prec
		if (USE_PREC == 1) {
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');

			//cluster.domains[d].K.MatAddInPlace( cluster.domains[d]._RegMat, 'N', -1.0);
			cluster.domains[d].K.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			cluster.domains[d]._RegMat.MatVecCOO(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N', -1.0);
			//cluster.domains[d]._RegMat.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N', 0, 0, -1.0);

			//cluster.domains[d].K.MatAddInPlace( cluster.domains[d]._RegMat, 'N',  1.0);

			//cluster.domains[d].Kplus.ImportMatrix_wo_Copy (cluster.domains[d].K);

			//cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
		}

		if (USE_PREC == 2) { // Weight scaling function
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster2[d], 'T');
			//cluster.x_prim_cluster2[d] = cluster.x_prim_cluster1[d];
		}
	}

	std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
	SEQ_VECTOR < double > y_out_tmp;
	for (eslocal d = 0; d < cluster.domains.size(); d++) {
		y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
		cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

		for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling
	}
	time_eval.timeEvents[0].end();


	time_eval.timeEvents[1].start();
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[1].end();


	time_eval.totalTime.end();

}

// *** END - Iteration solver class *************************************
// **********************************************************************



// **********************************************************************
// *** Communication layer **********************************************
void   SendMatrix  ( eslocal rank, eslocal source_rank, SparseMatrix & A_in, eslocal dest_rank, SparseMatrix & B_out) {

	eslocal param_tag = 1;
	eslocal I_row_tag = 2;
	eslocal J_col_tag = 3;
	eslocal V_val_tag = 4;

	MPI_Status status;
	MPI_Request request;

	if (rank == source_rank) {
		eslocal send_par_buf[4];
		send_par_buf[0] = A_in.cols;
		send_par_buf[1] = A_in.rows;
		send_par_buf[2] = A_in.nnz;
		send_par_buf[3] = A_in.type;

#ifdef XE6
		MPI_Send(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, MPI_COMM_WORLD);
		MPI_Send(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, MPI_COMM_WORLD );
		MPI_Send(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, MPI_COMM_WORLD );
		MPI_Send(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD );
#else
		MPI_Isend(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD, & request);
#endif

	}

	if (rank == dest_rank) {
		eslocal recv_par_buf[4];
		MPI_Recv(recv_par_buf, 4, esglobal_mpi, source_rank, param_tag, MPI_COMM_WORLD, & status);
		B_out.cols = recv_par_buf[0];
		B_out.rows = recv_par_buf[1];
		B_out.nnz  = recv_par_buf[2];
		B_out.type = recv_par_buf[3];

		B_out.CSR_I_row_indices.resize(B_out.rows + 1);
		B_out.CSR_J_col_indices.resize(B_out.nnz);
		B_out.CSR_V_values.     resize(B_out.nnz);

		MPI_Recv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, esglobal_mpi,    source_rank, I_row_tag, MPI_COMM_WORLD, & status );
		MPI_Recv(&B_out.CSR_J_col_indices[0], B_out.nnz,      esglobal_mpi,    source_rank, J_col_tag, MPI_COMM_WORLD, & status );
		MPI_Recv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE, source_rank, V_val_tag, MPI_COMM_WORLD, & status );
	}

#ifdef WIN32
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void   SendMatrix ( SparseMatrix & A_in, eslocal dest_rank ) {

	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);


	eslocal param_tag = 1;
	eslocal I_row_tag = 2;
	eslocal J_col_tag = 3;
	eslocal V_val_tag = 4;

	MPI_Status status;
	MPI_Request request;

	eslocal send_par_buf[4];
	send_par_buf[0] = A_in.cols;
	send_par_buf[1] = A_in.rows;
	send_par_buf[2] = A_in.nnz;
	send_par_buf[3] = A_in.type;

//#ifdef XE6
//		MPI_Send(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, MPI_COMM_WORLD);
//		MPI_Send(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, MPI_COMM_WORLD );
//		MPI_Send(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, MPI_COMM_WORLD );
//		MPI_Send(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD );
//#else
		MPI_Isend(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, MPI_COMM_WORLD, & request);
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, MPI_COMM_WORLD, & request);
//#endif

#ifdef WIN32
//	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void   RecvMatrix ( SparseMatrix & B_out, eslocal source_rank) {

	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);


	eslocal param_tag = 1;
	eslocal I_row_tag = 2;
	eslocal J_col_tag = 3;
	eslocal V_val_tag = 4;

	MPI_Status status;
	MPI_Request request;

	eslocal recv_par_buf[4];
	MPI_Recv(recv_par_buf, 4, esglobal_mpi, source_rank, param_tag, MPI_COMM_WORLD, & status);
	B_out.cols = recv_par_buf[0];
	B_out.rows = recv_par_buf[1];
	B_out.nnz  = recv_par_buf[2];
	B_out.type = recv_par_buf[3];

	B_out.CSR_I_row_indices.resize(B_out.rows + 1);
	B_out.CSR_J_col_indices.resize(B_out.nnz);
	B_out.CSR_V_values.     resize(B_out.nnz);

	MPI_Recv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, esglobal_mpi,    source_rank, I_row_tag, MPI_COMM_WORLD, & status );
	MPI_Recv(&B_out.CSR_J_col_indices[0], B_out.nnz,      esglobal_mpi,    source_rank, J_col_tag, MPI_COMM_WORLD, & status );
	MPI_Recv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE, source_rank, V_val_tag, MPI_COMM_WORLD, & status );


#ifdef WIN32
	//MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void   BcastMatrix ( eslocal rank, eslocal mpi_root, eslocal source_rank, SparseMatrix & A) {

	//eslocal param_tag = 1;
	//eslocal I_row_tag = 2;
	//eslocal J_col_tag = 3;
	//eslocal V_val_tag = 4;

	//MPI_Status status;
	//MPI_Request request;

	eslocal send_par_buf[4];
	//eslocal recv_par_buf[4];

	if (rank == source_rank) {
		send_par_buf[0] = A.cols; send_par_buf[1] = A.rows; send_par_buf[2] = A.nnz; send_par_buf[3] = A.type;
	}

	MPI_Bcast(send_par_buf, 4, esglobal_mpi, source_rank, MPI_COMM_WORLD);

	if (rank != source_rank) {
		A.cols = send_par_buf[0]; A.rows = send_par_buf[1]; A.nnz  = send_par_buf[2]; A.type = send_par_buf[3];
		A.CSR_I_row_indices.resize(A.rows + 1);
		A.CSR_J_col_indices.resize(A.nnz);
		A.CSR_V_values.     resize(A.nnz);
	}

	MPI_Bcast(&A.CSR_I_row_indices[0], A.rows + 1, esglobal_mpi, source_rank, MPI_COMM_WORLD);
	MPI_Bcast(&A.CSR_J_col_indices[0], A.nnz,      esglobal_mpi, source_rank, MPI_COMM_WORLD);
	MPI_Bcast(&A.CSR_V_values[0],      A.nnz,   MPI_DOUBLE, source_rank, MPI_COMM_WORLD);
}

void   All_Reduce_lambdas_compB( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out )
{
	for (eslocal i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (eslocal j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			cluster.my_comm_lambdas[i][j] = x_in[cluster.my_comm_lambdas_indices_comp[i][j]];
		}
	}


	MPI_Request * mpi_req  = new MPI_Request [cluster.my_neighs.size()];
	MPI_Status  * mpi_stat = new MPI_Status  [cluster.my_neighs.size()];

	cluster.iter_cnt_comm++;
	eslocal tag = cluster.iter_cnt_comm;

	for (eslocal neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		MPI_Sendrecv(
			&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,
			&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag, MPI_COMM_WORLD, &mpi_stat[neigh_i] );
	}

	//for (eslocal neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
	//	eslocal b_size = cluster.my_comm_lambdas[neigh_i].size();
	//	MPI_Isend(&b_size,                              1                                      , esglobal_mpi   , cluster.my_neighs[neigh_i], tag + 100, MPI_COMM_WORLD, &mpi_req[neigh_i] );
	//	MPI_Isend(&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,       MPI_COMM_WORLD, &mpi_req[neigh_i] );
	//
	//}

	//for (eslocal neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
	//	eslocal r_size = 0;
	//	MPI_Recv(&r_size                             ,                                       1, esglobal_mpi   , cluster.my_neighs[neigh_i], tag + 100, MPI_COMM_WORLD, &mpi_stat[neigh_i] );
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
	for (eslocal i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (eslocal j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			y_out[cluster.my_comm_lambdas_indices_comp[i][j]] += cluster.my_recv_lambdas[i][j];
		}
	}



}

void   compress_lambda_vector  ( Cluster & cluster, SEQ_VECTOR <double> & decompressed_vec_lambda)
{
	//compress vector for CG in main loop
	for (eslocal i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[cluster.my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(cluster.my_lamdas_indices.size());
}

void   decompress_lambda_vector( Cluster & cluster, SEQ_VECTOR <double> & compressed_vec_lambda)
{
	SEQ_VECTOR <double> decompressed_vec_lambda (cluster.domains[0].B1.rows,0);

	for (eslocal i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[cluster.my_lamdas_indices[i]] = compressed_vec_lambda[i];

	compressed_vec_lambda = decompressed_vec_lambda;
}

double parallel_norm_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector )
{

	double wl = 0; double wg = 0;

	for (eslocal i = 0; i < cluster.my_lamdas_indices.size(); i++)
		wl = wl + (input_vector[i] * input_vector[i] * cluster.my_lamdas_ddot_filter[i]);

	MPI_Allreduce( &wl, &wg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double norm_l = sqrt(wg);

	return norm_l;
}

double parallel_ddot_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector1, SEQ_VECTOR<double> & input_vector2 )
{
	double a1 = 0; double a1g = 0;

	for (eslocal i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
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

	for (eslocal i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
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
