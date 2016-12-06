//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#include "itersolver.h"

using namespace espreso;

IterSolverBase::IterSolverBase():
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


void IterSolverBase::Setup ( SEQ_VECTOR <double> & parameters , Cluster & cluster_in )
{

	// *** MPI variables  **********************************************************
	//mpi_rank;	//mpi_root;	//mpi_size;

	MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	/* get number of processes */
	mpi_root = 0;

}

void IterSolverBase::Preprocessing ( Cluster & cluster )
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

void IterSolverBase::Solve_singular ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel )
{

	switch (config::solver::CGSOLVER) {
	case config::solver::CGSOLVERalternative::STANDARD:
		Solve_RegCG_singular_dom ( cluster, in_right_hand_side_primal );
		break;
	case config::solver::CGSOLVERalternative::PIPELINED:
		Solve_PipeCG_singular_dom( cluster, in_right_hand_side_primal );
		break;
	case config::solver::CGSOLVERalternative::FULL_ORTOGONAL:
		Solve_full_ortho_CG_singular_dom (cluster, in_right_hand_side_primal );
		break;
	case config::solver::CGSOLVERalternative::GMRES:
		Solve_GMRES_singular_dom (cluster, in_right_hand_side_primal );
		break;
	case config::solver::CGSOLVERalternative::BICGSTAB:
		Solve_BICGSTAB_singular_dom(cluster, in_right_hand_side_primal );
		break;
	case config::solver::CGSOLVERalternative::QPCE:
		Solve_QPCE_singular_dom(cluster, in_right_hand_side_primal );
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown CG solver";
	}

	 postproc_timing.totalTime.start();

	 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
	 timeGetSol.start();
	GetSolution_Primal_singular_parallel( cluster, in_right_hand_side_primal, out_primal_solution_parallel );
	 timeGetSol.endWithBarrier();
	 postproc_timing.addEvent(timeGetSol);

	 postproc_timing.totalTime.endWithBarrier();




}

void IterSolverBase::Solve_non_singular ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel )
{
	switch (config::solver::CGSOLVER) {
	case config::solver::CGSOLVERalternative::STANDARD:
		Solve_RegCG_nonsingular  ( cluster, in_right_hand_side_primal, out_primal_solution_parallel);
		break;
	case config::solver::CGSOLVERalternative::PIPELINED:
		Solve_PipeCG_nonsingular ( cluster, in_right_hand_side_primal, out_primal_solution_parallel);
		break;
//	case config::solver::CGSOLVERalternative::FULL_ORTOGONAL:
//		break;
//	case config::solver::CGSOLVERalternative::GMRES:
//		break;
//	case config::solver::CGSOLVERalternative::BICGSTAB:
//		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown CG solver";
	}
}


void IterSolverBase::GetResiduum_Dual_singular_parallel    ( Cluster & cluster, SEQ_VECTOR <double> & dual_residuum_out ) {

	dual_residuum_out = dual_residuum_compressed_parallel;
	cluster.decompress_lambda_vector( dual_residuum_out );

}

void IterSolverBase::GetSolution_Dual_singular_parallel    ( Cluster & cluster, SEQ_VECTOR <double> & dual_solution_out, SEQ_VECTOR<double> & amplitudes_out ) {

	dual_solution_out = dual_soultion_compressed_parallel;
	cluster.decompress_lambda_vector( dual_solution_out );

	amplitudes_out	  = amplitudes;

}

void IterSolverBase::GetSolution_Primal_singular_parallel  ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
		SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out ) {

	MakeSolution_Primal_singular_parallel(cluster, in_right_hand_side_primal, primal_solution_out );
	//primal_solution_out = primal_solution_parallel;


	// KKT conditions
	// TODO: OPTIMIZATION

	bool check_solution = true;
	if (check_solution) {

		size_t dl_size = cluster.my_lamdas_indices.size();

		eslocal ieq_size = 0;
		eslocal eq_size = 0; //cluster.my_lamdas_indices.size();

		SEQ_VECTOR <double> KKT1_norm (cluster.domains.size(), 0.0);
		double KKT1_norm_cluster_local = 0.0;
		double KKT1_norm_cluster_global = 0.0;

		SEQ_VECTOR <double> KKT1_norm2 (cluster.domains.size(), 0.0);
		double KKT1_norm_cluster_local2 = 0.0;
		double KKT1_norm_cluster_global2 = 0.0;

		// KKT 1

		for (size_t d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
			for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = dual_soultion_compressed_parallel[ cluster.domains[d].lambda_map_sub_local[i]]; // * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T'); // Bt*lambda
		}

		for (size_t d = 0; d < cluster.domains.size(); d++) {
			cluster.domains[d].K.MatVec(primal_solution_out[d], cluster.x_prim_cluster2[d],'N');
			cluster.domains[d]._RegMat.MatVecCOO(primal_solution_out[d], cluster.x_prim_cluster2[d],'N', 1.0, -1.0); // K*u
		}

		for (size_t d = 0; d < cluster.domains.size(); d++) {
			for (size_t pi = 0; pi < primal_solution_out[d].size(); pi++ ) {
				cluster.x_prim_cluster1[d][pi] =   in_right_hand_side_primal[d][pi]
										         - cluster.x_prim_cluster1[d][pi];		// f - Bt*lambda


				KKT1_norm2[d] += cluster.x_prim_cluster1[d][pi] * cluster.x_prim_cluster1[d][pi]; // norm (f - Bt*lambda)

				cluster.x_prim_cluster1[d][pi] = cluster.x_prim_cluster2[d][pi] - cluster.x_prim_cluster1[d][pi]; // K*u - (f - bt*lambda)
				KKT1_norm[d] += cluster.x_prim_cluster1[d][pi] * cluster.x_prim_cluster1[d][pi]; // norm (K*u - (f - bt*lambda))
			}
			KKT1_norm_cluster_local += KKT1_norm[d];
			KKT1_norm_cluster_local2 += KKT1_norm2[d];

		}

		MPI_Allreduce( &KKT1_norm_cluster_local, &KKT1_norm_cluster_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		KKT1_norm_cluster_global = sqrt(KKT1_norm_cluster_global);

		MPI_Allreduce( &KKT1_norm_cluster_local2, &KKT1_norm_cluster_global2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		KKT1_norm_cluster_global2 = sqrt(KKT1_norm_cluster_global2);

		//ESINFO(CONVERGENCE) << " KKT1 norm:              " << std::setw(6) << KKT1_norm_cluster_global / KKT1_norm_cluster_global2;


		// KKT2

		SEQ_VECTOR <double> vec_c_l   (dl_size, 0);
		SEQ_VECTOR <double> lb        (dl_size, 0);
		SEQ_VECTOR <double> Bu_l      (dl_size, 0.0);
		SEQ_VECTOR <double> Be_l      (dl_size, 0.0);
		SEQ_VECTOR <double> Bn_l      (dl_size, 0.0);
		SEQ_VECTOR <double> Bn_lLambda (dl_size, 0.0);
		SEQ_VECTOR <double> lambdan_l (dl_size, 0.0);
		SEQ_VECTOR <double> ce_l      (dl_size, 0.0);
		SEQ_VECTOR <double> cn_l      (dl_size, 0.0);

		double max_Bn_l = 0.0;
		double max_Bn_l_g = 0.0;

		double norm_ce = 0.0;
		double norm_cn = 0.0;
		double norm_Beu = 0.0;
		double norm_Bn_lLambda = 0.0;

		double lambda_n_max = -INFINITY;
		double lambda_n_max_g = 0.0;

		double lambda_n_max_2 = 0.0;
		double lambda_n_max_2_g = 0.0;


		cluster.CreateVec_c_perCluster ( vec_c_l );
		cluster.CreateVec_lb_perCluster ( lb );

		for (size_t i = 0; i < cluster.compressed_tmp.size(); i++)
			cluster.compressed_tmp[i] = 0.0;

		for (size_t d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR <double> tmp_dual(cluster.domains[d].B1_comp_dom.rows, 0.0);
			cluster.domains[d].B1_comp_dom.MatVec (primal_solution_out[d], tmp_dual, 'N', 0, 0, 0.0);
			for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += tmp_dual[i];
				//Bu_l[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.x_prim_cluster1[d][i];
			}
		}

		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, Bu_l);

		for (size_t i = 0; i < vec_c_l.size(); i++){
			if ( lb[i] == 0 ) {
				lambdan_l[i] = dual_soultion_compressed_parallel[i];

				if (lambda_n_max < lambdan_l[i])
					lambda_n_max = lambdan_l[i];

				if (lambda_n_max_2 < -lambdan_l[i])
					lambda_n_max_2 = -lambdan_l[i];

				cn_l[i] = vec_c_l [i];
				Bn_l[i] = Bu_l[i] - vec_c_l [i];

				if (max_Bn_l < Bn_l[i])
					max_Bn_l = Bn_l[i];

				Bn_lLambda[i] = Bn_l[i] * lambdan_l[i];

				ieq_size++;
			} else {
				ce_l[i] = vec_c_l[i];
				Be_l[i] = Bu_l[i] - vec_c_l [i];
				eq_size++;
			}
		}

		MPI_Allreduce( &max_Bn_l, &max_Bn_l_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		max_Bn_l_g = fabs(max_Bn_l_g);

		MPI_Allreduce( &lambda_n_max, &lambda_n_max_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (fabs(lambda_n_max_g) < 10e-8)
			lambda_n_max_g += 1.0;


		MPI_Allreduce( &lambda_n_max_2, &lambda_n_max_2_g, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		lambda_n_max_2_g = fabs(lambda_n_max_2_g);

		norm_ce  = parallel_norm_compressed(cluster, ce_l);
		if (fabs(norm_ce) < 10e-8)
			norm_ce += 1.0;

		norm_cn  = parallel_norm_compressed(cluster, cn_l);
		if (fabs(norm_cn) < 10e-8)
			norm_cn += 1.0;

		norm_Beu = parallel_norm_compressed(cluster, Be_l);
		parallel_norm_compressed(cluster, Bn_l);
		norm_Bn_lLambda = parallel_norm_compressed(cluster, Bn_lLambda);



		ESINFO(CONVERGENCE) << " **** Karush-Kuhn-Tucker-conditions ****     ";


		ESINFO(CONVERGENCE) << " Solution norm: norm(K*u - f + Bt*Lambda) =      				" << std::setw(6) << KKT1_norm_cluster_global / KKT1_norm_cluster_global2;

		ESINFO(CONVERGENCE) << " Equality constraints: norm(Be*u - ce) =					" << std::setw(6) << norm_Beu / norm_ce;

		// KKT3
		ESINFO(CONVERGENCE) << " Inequality constraints: norm(max(Bn*u - cn,0)) =				" << std::setw(6) << max_Bn_l_g / norm_cn;
		ESINFO(CONVERGENCE) << " Check Multipliers positiveness: norm(max(-Lambda_N,0) =			" << std::setw(6) << lambda_n_max_2_g / lambda_n_max_g;
		ESINFO(CONVERGENCE) << " Check norm of Normal Multipliers: norm((Be*u - ce)*Lambda_N) =			" << std::setw(6) << norm_Bn_lLambda / ( norm_cn * lambda_n_max_g );

//		switch (config::solver::CGSOLVER) {
//		case config::solver::CGSOLVERalternative::QPCE:
//			Solve_QPCE_singular_dom(cluster, in_right_hand_side_primal );
//			break;
//		default:
//			ESINFO(GLOBAL_ERROR) << "Unknown CG solver";
//		}


	}
}

void IterSolverBase::MakeSolution_Primal_singular_parallel ( Cluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
		SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out )  {

	//primal_solution_parallel.clear();
	primal_solution_out.clear();

	// R * mu
	SEQ_VECTOR<SEQ_VECTOR<double> > R_mu_prim_cluster;

	for (size_t d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR <double > tmp (cluster.domains[d].domain_prim_size);
		if (USE_HFETI == 1)

			if ( config::solver::REGULARIZATION == config::solver::REGULARIZATIONalternative::FIX_POINTS )
				cluster.domains[d].Kplus_R.DenseMatVec(amplitudes, tmp, 'N', 0, 0);
		  	else
		  		cluster.domains[d].Kplus_Rb.DenseMatVec(amplitudes, tmp, 'N', 0, 0);

		else
			cluster.domains[d].Kplus_R.DenseMatVec(amplitudes, tmp, 'N', d * cluster.domains[d].Kplus_R.cols, 0);

		R_mu_prim_cluster.push_back(tmp);
	}

	for (size_t d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
		SEQ_VECTOR < double > tmp      ( cluster.domains[d].domain_prim_size  );

		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = dual_soultion_compressed_parallel[ cluster.domains[d].lambda_map_sub_local[i]];

		cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, tmp, 'T');

		for (size_t i = 0; i < tmp.size(); i++)
			tmp[i] = in_right_hand_side_primal[d][i] - tmp[i];
			//tmp[i] = cluster.domains[d].f[i] - tmp[i];

		//primal_solution_parallel.push_back(tmp);
		primal_solution_out.push_back(tmp);
	}

	if ( cluster.USE_HFETI == 0) {
		for (size_t d = 0; d < cluster.domains.size(); d++)
			cluster.domains[d].multKplusLocal(primal_solution_out[d]);
			//cluster.domains[d].multKplusLocal(primal_solution_parallel[d]);
	} else  {
		//cluster.multKplusGlobal_l(primal_solution_parallel);
		cluster.multKplusGlobal_l(primal_solution_out);
	}

 	for (size_t d = 0; d < cluster.domains.size(); d++) {
		//for (size_t i = 0; i < primal_solution_parallel[d].size()	; i++) {
 		for (size_t i = 0; i < primal_solution_out[d].size()	; i++) {
 			primal_solution_out[d][i] = primal_solution_out[d][i] + R_mu_prim_cluster[d][i];
			//primal_solution_parallel[d][i] = primal_solution_parallel[d][i] + R_mu_prim_cluster[d][i];
			//primal_solution_parallel[d][i] = cluster.domains[d].up0[i] + R_mu_prim_cluster[d][i];
		}
	}

}


// POWER Method
double IterSolverBase::Solve_power_method ( Cluster & cluster, double tol, eslocal maxit, eslocal method)
{
	size_t dl_size = cluster.my_lamdas_indices.size();
	double norm_V_0 = 0;
	double err = 1;
	double lambda = 0;
	double lambda0 = 0;
	eslocal nit = 0;
	SEQ_VECTOR <double> V_0 (dl_size, 0);
	SEQ_VECTOR <double> V (dl_size, 0);
	SEQ_VECTOR <double> Y (dl_size, 0);
	SEQ_VECTOR <double> X (dl_size, 0);
	SEQ_VECTOR <double> Z (dl_size, 0);

    // 1 -1 1 -1 1 -1 ....... local global mapping
	for ( size_t i=0; i< dl_size; i++ )
	{
		if ( i%2 == 0 )
		{
			V_0[i] = 1;
		}
		else
		{
			V_0[i] = -1;
		}
	}

	norm_V_0 = parallel_norm_compressed(cluster, V_0);

	for ( size_t i=0; i< dl_size; i++ )
	{
		Y[i] = V_0[i]/norm_V_0;
	}



	if (method == 0){

	if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, Y, V , 0);
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, Y, V , 0);
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Y);


		if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, Y, V, 0);
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, Y, V, 0);
		}

	}else{

					if (USE_GGtINV == 1) {
						Projector_l_inv_compG( timeEvalProj, cluster, Y, V , 0);
					} else {
						Projector_l_compG    ( timeEvalProj, cluster, Y, V , 0);
					}

					apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Z);

					if (USE_GGtINV == 1) {
						Projector_l_inv_compG( timeEvalProj, cluster, Z, X, 0);
					} else {
						Projector_l_compG    ( timeEvalProj, cluster, Z, X, 0);
					}

					for (size_t i = 0; i < Z.size(); i++){
						V[i] = X[i]/method + ( Y[i]-V[i]);
					}
	    		}


//	apply_A_l_comp_dom_B(timeEvalAppa, cluster, Y, V);
    lambda = parallel_norm_compressed(cluster, V);

    for ( eslocal i=0; i < maxit; i++)
    {
    	if (err < tol )
    	{
    		break;
    	}
    	for ( size_t i=0; i< dl_size; i++ )
    	{
    		Y[i] = V[i]/lambda;
    	}


    	if (method == 0){

    		if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, Y, V , 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, Y, V , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Y);


				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, Y, V, 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, Y, V, 0);
				}

    	}else{

				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, Y, V , 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, Y, V , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Z);

				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, Z, X, 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, Z, X, 0);
				}

				for (size_t i = 0; i < Z.size(); i++){
					V[i] = X[i]/method + ( Y[i]-V[i]);
				}
    		}

    	//apply_A_l_comp_dom_B(timeEvalAppa, cluster, Y, V);
    	lambda0 = lambda;
    	lambda = parallel_norm_compressed(cluster, V);
    	err = fabs(lambda - lambda0)/fabs(lambda);
    	nit++;

    }
    // TODO return number of iterations
    return lambda;
}

//	static bool abs_compare(int a, int b)
//	{
//	    return (std::abs(a) < std::abs(b));
//	}

void IterSolverBase::proj_gradient ( SEQ_VECTOR <double> & x,
		SEQ_VECTOR <double> & g,
		SEQ_VECTOR <double> & lb,
		double alpha, double prec, SEQ_VECTOR <double> & g_til, SEQ_VECTOR <double> & fi_til, SEQ_VECTOR <double> & beta_til,
		SEQ_VECTOR <bool> & free )
{

//	double norm_x_l =0;
//	double norm_x =0;
//
//	std::vector<eslocal>::iterator result;
//	norm_x_l = *std::max_element(x.begin(), x.end(), abs_compare);
//
//	norm_x_l = fabs(norm_x_l);
//
//	MPI_Allreduce(&norm_x_l, &norm_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	for ( size_t i=0; i< x.size(); i++ )
	{
		g_til[i] = 1.0/alpha * (x[i]-( std::max(x[i] - alpha * g[i], lb[i] ) ) );

		if ( g_til[i] * g[i] < 0 )
		{
			g_til[i] = 0;
		}

		free[i] = ( x[i] - lb[i] ) > prec;//(prec * norm_x);


		fi_til[i] = free[i] * g_til[i];
		beta_til[i] = !free[i] * g_til[i];
	}

}




// QPCE is a variant of the algorithm SMALBE   (SemiMonotonous Augmented Lagrangians with Bound and Equality constraints)
void IterSolverBase::Solve_QPCE_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{

	double _epsilon = 1e-4;
	eslocal _maxit = 100;
	eslocal _maxit_in = 200;
	double _Gamma = 1;
	double _M = 1;
	double _rho = 1;
	double _eta = 1;
	double _beta = 0.8;
	double _alpham = 2;
	double _precQ = 1e-12;
	double _epsilon_power = 1e-8;
	eslocal _maxit_power = 55;
	eslocal _method = 0;
	eslocal halfStep = 0;

	eslocal output_n_it = 0;
	eslocal output_n_it_in = 0;
	eslocal output_n_cg = 0;
	eslocal output_n_prop = 0;
	eslocal output_n_exp = 0;
	eslocal output_n_hess = 0;

	eslocal sum_output_n_cg = 0;
	eslocal sum_output_n_prop = 0;
	eslocal sum_output_n_exp = 0;
	eslocal sum_output_n_hess = 0;

	eslocal dl_size = cluster.my_lamdas_indices.size();
	SEQ_VECTOR <double> b_l  (dl_size, 0);
	SEQ_VECTOR <double> b_l_  (dl_size, 0);
	SEQ_VECTOR <double> x_im  (dl_size, 0);
	SEQ_VECTOR <double> Ax_im  (dl_size, 0);
	SEQ_VECTOR <double> tmp  (dl_size, 0);
	SEQ_VECTOR <double> lb (dl_size, 0);
	SEQ_VECTOR <double> x_l (dl_size, 0);
	SEQ_VECTOR <double> bCtmu (dl_size, 0);
	SEQ_VECTOR <double> g_l (dl_size, 0);
	SEQ_VECTOR <double> Ax_l (dl_size, 0);
	SEQ_VECTOR <double> PAPx_l (dl_size, 0);
	SEQ_VECTOR <double> r_l (dl_size, 0);


	SEQ_VECTOR <double> w_l (dl_size, 0);
	SEQ_VECTOR <double> y_l (dl_size, 0);
	SEQ_VECTOR <double> tmp_2 (dl_size, 0);



	SEQ_VECTOR <double> nx_l (dl_size, 0);
	SEQ_VECTOR <double> ng_l (dl_size, 0);



	SEQ_VECTOR <double> mu (cluster.G1_comp.rows, 0);
	SEQ_VECTOR <double> mu_tmp (cluster.G1_comp.rows, 0);


	SEQ_VECTOR <double> g_til (dl_size, 0);
	SEQ_VECTOR <double> fi_til (dl_size, 0);
	SEQ_VECTOR <double> beta_til (dl_size, 0);
	SEQ_VECTOR <bool> _free (dl_size, 0);
	SEQ_VECTOR <double> p_l (dl_size, 0);
	SEQ_VECTOR <double> Cx_l (cluster.G1_comp.rows, 0);
	double normx_l = 0;
	SEQ_VECTOR <double> test_vec (dl_size, 0);

	double beta_til_g_l = 0;
	double fi_til_g_l = 0;

	double pAp = 0;
	double pg = 0;
	double alpha_cg = 0;
	double alpha_f_l = 0;
	double alpha_f = 0;
	eslocal cnt_l = 0;
	eslocal cnt = 0;
	double gamma_p = 0;
	double norm_test_vec = 0;
	double normCx = 0;
	double normCx_x = 0;
	SEQ_VECTOR <double> bCtmu_prev (dl_size, 0);

	double lag0 = -INFINITY;
	double lag1 = 0;

	double maxeig = Solve_power_method ( cluster, _epsilon_power, _maxit_power, _method);

//	double maxeig_old = 0.0;
//
//	for (eslocal i = 0; i < 10; i++){
//
//		maxeig = Solve_power_method ( cluster, _epsilon_power, _maxit_power, maxeig_old);
//
//		if ( std::fabs( maxeig-maxeig_old ) < 1e-8){
//
//			break;
//
//		} else {
//			maxeig_old = maxeig;
//		}
//
//	}

//	double alpha = _alpham/maxeig;
//	double rho = _rho*maxeig;
//	maxeig = 1.0;

	double alpha = _alpham;
	double rho = _rho;

	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	cluster.CreateVec_lb_perCluster ( lb );

	// BEGIN*** projection of right hand side b
	for (size_t i = 0; i < b_l.size(); i++){
		b_l_[i] = b_l[i];
	}

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, b_l_, b_l , 0);
	} else {
		Projector_l_compG    ( timeEvalProj, cluster, b_l_, b_l , 0);
	}
	// END*** projection of right hand side b



	// BEGIN*** Homogenization of the equality constraints and initialization
	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_im, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_im, 1 );
	}

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_im, Ax_im);

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, Ax_im, tmp , 0);
	} else {
		Projector_l_compG    ( timeEvalProj, cluster, Ax_im, tmp , 0);
	}

	for (size_t i = 0; i < tmp.size(); i++){
		b_l[i] = b_l[i] - tmp[i];
		lb[i] = lb[i] - x_im[i];
		x_l[i] = std::max( lb[i] , 0.0 );

		b_l[i] = b_l[i]/maxeig;


		bCtmu[i] = b_l[i];
	}


	//double norm_b = parallel_norm_compressed(cluster, b_l);
	//double tol = _epsilon * norm_b;
	// END*** Homogenization of the equality constraints and initialization

    /// BEGIN*** Hessian PAP+rho*Ct*inv(C*Ct)*C
//	if (USE_GGtINV == 1) {
//		Projector_l_inv_compG( timeEvalProj, cluster, x_l, tmp , 0);
//	} else {
//		Projector_l_compG    ( timeEvalProj, cluster, x_l, tmp , 0);
//	}
//
//	apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);
//
//	for (size_t i = 0; i < tmp.size(); i++){
//		Ax_l[i] = Ax_l[i] - rho * x_l[i];
//	}
//
//	if (USE_GGtINV == 1) {
//		Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
//	} else {
//		Projector_l_compG    ( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
//	}
//
//	for (size_t i = 0; i < tmp.size(); i++){
//		PAPx_l[i] = PAPx_l[i] + rho * x_l[i];
//	}



		if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, x_l, tmp , 0);
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, x_l, tmp , 0);
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

		if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
		}

		for (size_t i = 0; i < tmp.size(); i++){
			PAPx_l[i] = PAPx_l[i]/maxeig + rho * ( x_l[i]-tmp[i]);
		}


	// END*** Hessian PAP+rho*Ct*inv(C*Ct)*C
	sum_output_n_hess++;

//	for (eslocal i = 0; i < x_l.size(); i++){
//		x_l[i] = std::max( lb[i], 0.0 );
//	}


	for (size_t i = 0; i < tmp.size(); i++){
		g_l[i] = PAPx_l[i] - bCtmu[i];
	}

	proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );

	for (size_t i = 0; i < tmp.size(); i++){
		p_l[i] = _free[i] * g_l[i];
		//test_vec[i] = g_til[i];
	}







	switch (USE_PREC) {
			case config::solver::PRECONDITIONERalternative::LUMPED:
			case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
			case config::solver::PRECONDITIONERalternative::DIRICHLET:
			case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
			case config::solver::PRECONDITIONERalternative::MAGIC:

				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, p_l, w_l, 0 );
				} else {
					Projector_l_compG( timeEvalProj, cluster, p_l, w_l, 0 );
				}

				apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, tmp_2);

				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, tmp_2, y_l, 0 );
				} else {
					Projector_l_compG		  ( timeEvalProj, cluster, tmp_2, y_l, 0 );
				}

				for (size_t k = 0; k < p_l.size(); k++){
					p_l[k] = (y_l[k] * maxeig + 1.0 / rho * ( p_l[k]-w_l[k])) * _free[k];

				}

				break;
			case config::solver::PRECONDITIONERalternative::NONE:


				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
			}








	cluster.G1_comp.MatVec(x_l, Cx_l, 'N');

	normx_l = parallel_norm_compressed(cluster, x_l);

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, Cx_l, tmp, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, Cx_l, tmp, 1 );
	}

	normCx = sqrt( parallel_ddot_compressed(cluster, x_l, tmp) );

	if (normx_l == 0){
		normCx_x = INFINITY;
	} else{
		normCx_x = normCx/normx_l;
	}

	double norm_b = parallel_norm_compressed(cluster, b_l);
	double tol = _epsilon * norm_b;
	norm_test_vec = parallel_norm_compressed(cluster, g_til);

	double mchange = 0.0;


	ESINFO(CONVERGENCE) << "===================================================================================================";
	ESINFO(CONVERGENCE) << "	QUADRATIC PROGRAMMING WITH SIMPLE BOUNDS AND EQUALITY CONSTRAINTS (QPCE)";
	ESINFO(CONVERGENCE) << "===================================================================================================";
	//ESINFO(CONVERGENCE) << "'Variables/equality constraints: n/m = %d/%d\n',n,m);
	ESINFO(CONVERGENCE) << "	Terminating tolerance: epsilon = "<< tol*maxeig ;
	//ESINFO(CONVERGENCE) << "Parameter settings:\n'); disp(options);
	ESINFO(CONVERGENCE) << "---------------------------------------------------------------------------------------------------";
	ESINFO(CONVERGENCE) << "Out_it   L(x,mu,rho)   ||~g(x)||       ||Cx||       Exp  Prop  Cgm   No_A      rho              M";
	ESINFO(CONVERGENCE) << "---------------------------------------------------------------------------------------------------";

	ESINFO(CONVERGENCE)<< std::setw(3) << output_n_it << std::setw(16) << lag1 << std::setw(15) << norm_test_vec*maxeig << std::setw(15)  <<normCx_x << std::setw(6)  <<  output_n_exp <<  std::setw(6) <<  output_n_prop <<  std::setw(6)  <<  output_n_cg <<  std::setw(6)  <<  output_n_hess <<  std::setw(15) <<  rho << std::setw(10) <<  _M ;

	eslocal stop = 0;
	//for (eslocal i=0; i < _maxit; i++){
	while ( stop == 0 ){
		output_n_it_in = 0;
		output_n_cg = 0;
		output_n_exp = 0;
		output_n_prop= 0;
		output_n_hess= 0;


		while ( (norm_test_vec > (std::min(_M * normCx ,_eta * norm_b)))  && (output_n_it_in < _maxit_in ) && (!((norm_test_vec <= tol) && (normCx <= _epsilon * normx_l))) ) {

			beta_til_g_l = parallel_ddot_compressed(cluster, beta_til, g_l);
			fi_til_g_l = parallel_ddot_compressed(cluster, fi_til, g_l);
			if ( std::max(0.0,beta_til_g_l) <=  _Gamma * std::max(0.0,fi_til_g_l) ){

				/// HESSS
//				if (USE_GGtINV == 1) {
//					Projector_l_inv_compG( timeEvalProj, cluster, p_l, tmp , 0);
//				} else {
//					Projector_l_compG    ( timeEvalProj, cluster, p_l, tmp , 0);
//				}
//
//				apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);
//
//				for (size_t k = 0; k < tmp.size(); k++) {
//					Ax_l[k] = Ax_l[k] - rho * p_l[k];
//				}
//
//				if (USE_GGtINV == 1) {
//					Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
//				} else {
//					Projector_l_compG    ( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
//				}
//
//				for (size_t k = 0; k < tmp.size(); k++) {
//					PAPx_l[k] = PAPx_l[k] + rho * p_l[k];
//				}

				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, p_l, tmp , 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, p_l, tmp , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
				}

				for (size_t k = 0; k < tmp.size(); k++){
					PAPx_l[k] = PAPx_l[k]/maxeig + rho * ( p_l[k]-tmp[k]);
				}
				/// HESSS
				output_n_hess++;
				sum_output_n_hess++;

				pAp = parallel_ddot_compressed(cluster, PAPx_l, p_l);
				pg = parallel_ddot_compressed(cluster, g_l, p_l);
				alpha_cg = pg/pAp;

				cnt_l = 0;
				for (size_t k=0; k < tmp.size(); k++) {
				   if ( p_l[k] > 0) {
					   tmp[cnt_l] = (x_l[k]-lb[k])/p_l[k];
					   cnt_l += 1;
				   }
				}

				MPI_Allreduce(&cnt_l, &cnt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

				if (cnt == 0) {
					alpha_f = INFINITY;
				} else {
					alpha_f_l = *std::min_element( tmp.begin(), tmp.begin() + cnt_l );
					MPI_Allreduce(&alpha_f_l, &alpha_f, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
				}

				if (alpha_cg <= alpha_f){

					for ( size_t k = 0; k < x_l.size(); k++){
						x_l[k] = x_l[k] - alpha_cg * p_l[k];
						g_l[k] = g_l[k] - alpha_cg * PAPx_l[k];
					}

					proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );

					for (size_t k = 0; k < tmp.size(); k++){
						tmp[k] = _free[k] * g_l[k];
					}



					switch (USE_PREC) {
							case config::solver::PRECONDITIONERalternative::LUMPED:
							case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
							case config::solver::PRECONDITIONERalternative::DIRICHLET:
							case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
							case config::solver::PRECONDITIONERalternative::MAGIC:

								if (USE_GGtINV == 1) {
									Projector_l_inv_compG( timeEvalProj, cluster, tmp, w_l, 0 );
								} else {
									Projector_l_compG( timeEvalProj, cluster, tmp, w_l, 0 );
								}

								apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, tmp_2);

								if (USE_GGtINV == 1) {
									Projector_l_inv_compG( timeEvalProj, cluster, tmp_2, y_l, 0 );
								} else {
									Projector_l_compG		  ( timeEvalProj, cluster, tmp_2, y_l, 0 );
								}

								for (size_t k = 0; k < tmp.size(); k++){
									tmp[k] = (y_l[k] * maxeig + 1.0 / rho * ( tmp[k]-w_l[k])) * _free[k];

								}



								break;
							case config::solver::PRECONDITIONERalternative::NONE:


								break;
							default:
								ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
							}





					gamma_p = parallel_ddot_compressed(cluster, tmp, PAPx_l)/pAp;

					for (size_t k = 0; k < p_l.size(); k++){
						p_l[k] = tmp[k] - gamma_p * p_l[k];
					}
					output_n_cg++;
					sum_output_n_cg++;

				} else {


					if ( halfStep == 1){

//						for ( size_t k = 0; k < nx_l.size(); k++){
//						//	nx_l[k] = std::max( lb[k], x_l[k] - alpha_cg * p_l[k]);
//							nx_l[k] = std::max( lb[k], x_l[k] - 10 * alpha_f * p_l[k]);
//						}




							for ( size_t k = 0; k < x_l.size(); k++){
								x_l[k] = x_l[k] -  (alpha_f) * p_l[k];
								g_l[k] = g_l[k] -  (alpha_f) * PAPx_l[k];
							}



						proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );

						for ( size_t k = 0; k < x_l.size(); k++){
							x_l[k] = x_l[k] - (alpha) * g_til[k];
						}

						if (USE_GGtINV == 1) {
							Projector_l_inv_compG( timeEvalProj, cluster, x_l, tmp , 0);
						} else {
							Projector_l_compG    ( timeEvalProj, cluster, x_l, tmp , 0);
						}

						apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

						if (USE_GGtINV == 1) {
							Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, g_l, 0);
						} else {
							Projector_l_compG    ( timeEvalProj, cluster, Ax_l, g_l, 0);
						}

						for (size_t k = 0; k < tmp.size(); k++){
							g_l[k] = g_l[k]/maxeig + rho * ( x_l[k]-tmp[k]);
						}
						/// HESSS
						output_n_hess++;
						sum_output_n_hess++;
						for (size_t k = 0; k < tmp.size(); k++){
							g_l[k] = g_l[k] - bCtmu[k];
						}

						proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );
						for (size_t k = 0; k < p_l.size(); k++){
							p_l[k] = _free[k] * g_l[k];
						}

						output_n_exp++;
						sum_output_n_exp++;



					} else {

						//alpha_f =alpha_f * 10;

						for ( size_t k = 0; k < x_l.size(); k++){
							x_l[k] = x_l[k] - alpha_f * p_l[k];
							g_l[k] = g_l[k] - alpha_f * PAPx_l[k];
						}

						proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );

						for ( size_t k = 0; k < x_l.size(); k++){
							x_l[k] = std::max( lb[k], x_l[k] - alpha_cg * g_til[k]);
						}

						if (USE_GGtINV == 1) {
							Projector_l_inv_compG( timeEvalProj, cluster, x_l, tmp , 0);
						} else {
							Projector_l_compG    ( timeEvalProj, cluster, x_l, tmp , 0);
						}

						apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

						if (USE_GGtINV == 1) {
							Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, g_l, 0);
						} else {
							Projector_l_compG    ( timeEvalProj, cluster, Ax_l, g_l, 0);
						}

						for (size_t k = 0; k < tmp.size(); k++){
							g_l[k] = g_l[k]/maxeig + rho * ( x_l[k]-tmp[k]);
						}
						/// HESSS
						output_n_hess++;
						sum_output_n_hess++;
						for (size_t k = 0; k < tmp.size(); k++){
							g_l[k] = g_l[k] - bCtmu[k];
						}

						proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );
						for (size_t k = 0; k < p_l.size(); k++){
							p_l[k] = _free[k] * g_l[k];
						}

						output_n_exp++;
						sum_output_n_exp++;
					}
				}
			} else {
				for (size_t k = 0; k < x_l.size(); k++){
					x_l[k] = x_l[k] -  alpha * g_til[k];
				}
				/// HESSS
//				if (USE_GGtINV == 1) {
//					Projector_l_inv_compG( timeEvalProj, cluster, x_l, tmp , 0);
//				} else {
//					Projector_l_compG    ( timeEvalProj, cluster, x_l, tmp , 0);
//				}
//
//				apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);
//
//				for (size_t k = 0; k < tmp.size(); k++) {
//					Ax_l[k] = Ax_l[k] - rho * x_l[k];
//				}
//
//				if (USE_GGtINV == 1) {
//					Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, g_l, 0);
//				} else {
//					Projector_l_compG    ( timeEvalProj, cluster, Ax_l, g_l, 0);
//				}
//
//				for (size_t k = 0; k < tmp.size(); k++) {
//					g_l[k] = g_l[k] + rho * x_l[k];
//				}
				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, x_l, tmp , 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, x_l, tmp , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

				if (USE_GGtINV == 1) {
					Projector_l_inv_compG( timeEvalProj, cluster, Ax_l, g_l, 0);
				} else {
					Projector_l_compG    ( timeEvalProj, cluster, Ax_l, g_l, 0);
				}

				for (size_t k = 0; k < tmp.size(); k++){
					g_l[k] = g_l[k]/maxeig + rho * ( x_l[k]-tmp[k]);
				}
				/// HESSS
				output_n_hess++;
				sum_output_n_hess++;

				for (size_t k = 0; k < tmp.size(); k++){
					g_l[k] = g_l[k] - bCtmu[k];
				}
				proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );

				for (size_t k = 0; k < tmp.size(); k++){
					p_l[k] = _free[k] * g_l[k];
				}

				output_n_prop++;
				sum_output_n_prop++;
			}

			output_n_it_in++;

			normx_l = parallel_norm_compressed(cluster, x_l);

			cluster.G1_comp.MatVec(x_l, Cx_l, 'N');

			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, Cx_l, tmp, 1 );
			} else {
				Projector_l_compG	 ( timeEvalProj, cluster, Cx_l, tmp, 1 );
			}

			normCx = sqrt( parallel_ddot_compressed(cluster, x_l, tmp) );

			//for (size_t k = 0; k < tmp.size(); k++){
			//	test_vec[k] = g_til[k];
			//}
			//norm_test_vec = parallel_norm_compressed(cluster, test_vec);

			norm_test_vec = parallel_norm_compressed(cluster, g_til);

		}
		output_n_it++;



		for (size_t k = 0; k < tmp.size(); k++) {
			tmp[k] = g_l[k] -bCtmu[k];
		}

		lag1  = parallel_ddot_compressed(cluster, tmp, x_l)*0.5;

		if (normx_l == 0){
			normCx_x = INFINITY;
		} else {
			normCx_x = normCx/normx_l;
		}

		ESINFO(CONVERGENCE)<< std::setw(3) << output_n_it << std::setw(16) << lag1 << std::setw(15) << norm_test_vec*maxeig << std::setw(15)  <<normCx_x << std::setw(6)  <<  output_n_exp <<  std::setw(6) <<  output_n_prop <<  std::setw(6)  <<  output_n_cg <<  std::setw(6)  <<  output_n_hess <<  std::setw(15) <<  rho << std::setw(10) <<  _M ;


		for (size_t k = 0; k < tmp.size(); k++) {
			bCtmu_prev[k] = bCtmu[k];
		}

		if ( ((norm_test_vec <= tol)  && (normCx <= _epsilon * normx_l)) || ( _maxit == output_n_it )) {
			stop = 1;
		} else {

			for (size_t k = 0; k < mu.size(); k++) {
				mu[k] = mu[k] + rho * Cx_l[k];
			}

			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, mu, tmp, 1 );
			} else {
				Projector_l_compG	 ( timeEvalProj, cluster, mu, tmp, 1 );
			}

			for (size_t k = 0; k < bCtmu.size(); k++) {
				bCtmu[k] = b_l[k] - tmp[k];
			}

			if ( !mchange && ( lag1 <= ( lag0 + rho* normCx*normCx*0.5 ))) {
				_M = _beta * _M; mchange = 1.0;
			} else {
				mchange = 0.0;
			}
			lag0 = lag1;
		}

		for (size_t k = 0; k < tmp.size(); k++) {
			g_l[k] = g_l[k] + bCtmu_prev[k]-bCtmu[k];
		}

		proj_gradient( x_l, g_l, lb, alpha, _precQ, g_til, fi_til, beta_til, _free );

		for (size_t k = 0; k < tmp.size(); k++) {
			p_l[k] = _free[k] * g_l[k];
		//	test_vec[k] = g_til[k];
		}
		//norm_test_vec = parallel_norm_compressed(cluster, test_vec);
		norm_test_vec = parallel_norm_compressed(cluster, g_til);
	}


	for (size_t k = 0; k < x_l.size(); k++) {
		x_l[k] = x_l[k] + x_im[k];
	}

	for (size_t k = 0; k < mu.size(); k++) {
		mu[k] = mu[k] * maxeig;
	}



	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, r_l);


	for (size_t k = 0; k < r_l.size(); k++) {
			r_l[k] = r_l[k] - b_l_[k];
		}

	cluster.G1_comp.MatVec(r_l, mu_tmp, 'N');

	for (size_t k = 0; k < mu.size(); k++) {
		mu_tmp[k] = mu_tmp[k] - mu[k];
		}

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, mu_tmp, amplitudes, 3 );
	} else {
		//TODO: Neni implementovan parametr 3
		Projector_l_compG	  ( timeEvalProj, cluster, mu_tmp, amplitudes, 3 );
	}

//
//	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, r_l);
//
//	for (eslocal k = 0; k < r_l.size(); k++) {
//		r_l[k] = r_l[k] - b_l_[k]-tmp[k];
//	}
//
	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;
//
//	if (USE_GGtINV == 1) {
//		Projector_l_inv_compG ( timeEvalProj, cluster, r_l, amplitudes, 2 );
//	} else {
//		Projector_l_compG	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
//	}
//
	for (size_t k = 0; k < amplitudes.size(); k++) {
		amplitudes[k] = -amplitudes[k];
	}


//	cnt_l = 0;	cnt = 0;	dnt_l = 0; 	dnt = 0;
//	for (eslocal k = 0; k < _free.size(); k++) {
//		cnt_l = cnt_l + _free[k];
//		dnt_l = dnt_l + !_free[k];
//	}
//	MPI_Allreduce(&cnt_l, &cnt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//	MPI_Allreduce(&dnt_l, &dnt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


	ESINFO(CONVERGENCE) << "---------------------------------------------------------------------------------------------------";
	//ESINFO(CONVERGENCE) << " Active/Free: " <<  std::setw(12) << cnt <<"/" << std::setw(12)  <<  dnt;
	ESINFO(CONVERGENCE) << " Expansion steps:            " << std::setw(6) << sum_output_n_exp;
	ESINFO(CONVERGENCE) << " Proportioning steps:        " << std::setw(6) << sum_output_n_prop;
	ESINFO(CONVERGENCE) << " CG steps:                   " << std::setw(6) << sum_output_n_cg;
//	ESINFO(CONVERGENCE) << " Balancing ratio updates:  %d\n', round(log(M/options.M)/log(options.beta)))";
	ESINFO(CONVERGENCE) << " Multiplications by Hessian: " << std::setw(6) << sum_output_n_hess;
	ESINFO(CONVERGENCE) << "===================================================================================================";
	ESINFO(CONVERGENCE) << "	END QPCE";
	ESINFO(CONVERGENCE) << "===================================================================================================";

}



void IterSolverBase::Solve_RegCG_singular_dom ( Cluster & cluster,
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
	for (size_t d = 0; d < cluster.domains.size(); d++)
		norm_prim_fl += cluster.domains[d].norm_f;

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Reduce   (&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	norm_prim_fg = sqrt(norm_prim_fg);



	// *** r = b - Ax *************************************************************
	#pragma omp parallel for
	for (size_t i = 0; i < r_l.size(); i++)
		r_l[i] = b_l[i] - Ax_l[i];

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, r_l, u_l , 0);
	} else {
		Projector_l_compG    ( timeEvalProj, cluster, r_l, u_l , 0);
	}

	// *** Calculate the stop condition *******************************************
	tol = epsilon * parallel_norm_compressed(cluster, u_l);

	int precision = ceil(log(1 / epsilon) / log(10)) + 1;
	int iterationWidth = ceil(log(CG_max_iter) / log(10));
	std::string indent = "   ";

	auto spaces = [] (int count) {
		std::stringstream ss;
		for (int i = 0; i < count; i++) {
			ss << " ";
		}
		return ss.str();
	};

	ESINFO(CONVERGENCE)
		<< spaces(indent.size() + iterationWidth - 4) << "iter"
		<< spaces(indent.size() + precision - 3) << "|r|" << spaces(2)
		<< spaces(indent.size() + 4) << "r" << spaces(4)
		<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "e" << spaces(precision / 2)
		<< spaces(indent.size()) << "time[s]";

	// *** Start the CG iteration loop ********************************************
	for (int iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.start();

		#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		switch (USE_PREC) {
		case config::solver::PRECONDITIONERalternative::LUMPED:
		case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
		case config::solver::PRECONDITIONERalternative::DIRICHLET:
		case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
		case config::solver::PRECONDITIONERalternative::MAGIC:
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
			break;
		case config::solver::PRECONDITIONERalternative::NONE:
			proj_time.start();
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector_l_compG		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj_time.end();

			#pragma omp parallel for
for (size_t i = 0; i < w_l.size(); i++)
				y_l[i] = w_l[i];

			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
		}


		//------------------------------------------
		if (iter == 0) {									// if outputs.n_it==1;

			#pragma omp parallel for
for (size_t i = 0; i < y_l.size(); i++)
				p_l[i] = y_l[i];							// p = y;

		} else {

			ddot_beta.start();
			beta_l =          parallel_ddot_compressed(cluster, y_l, w_l);
			beta_l = beta_l / parallel_ddot_compressed(cluster, yp_l, wp_l);
			ddot_beta.end();

			#pragma omp parallel for
for (size_t i = 0; i < p_l.size(); i++)
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
		#pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * Ap_l[i];
		}

		 norm_time.start();
		norm_l = parallel_norm_compressed(cluster, w_l);
		 norm_time.end();

		 timing.totalTime.end();

		ESINFO(CONVERGENCE)
			<< indent << std::setw(iterationWidth) << iter + 1
			<< indent << std::fixed << std::setprecision(precision) <<  norm_l / tol * epsilon
			<< indent << std::scientific << std::setprecision(3) << norm_l
			<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
			<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();

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

	switch (USE_PREC) {
	case config::solver::PRECONDITIONERalternative::LUMPED:
	case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
	case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
	case config::solver::PRECONDITIONERalternative::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case config::solver::PRECONDITIONERalternative::NONE:
		timing.addEvent(proj_time);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

	timing.addEvent(appA_time );
	timing.addEvent(ddot_beta);
	timing.addEvent(ddot_alpha);

	// *** END - Preslocal out the timing for the iteration loop ***********************************

}

void IterSolverBase::Solve_new_CG_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                            C G   -   N E W    I M P L E M E N T A T I O N                          # 
//##################################################################################################*/
//
	size_t dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);

	SEQ_VECTOR <double> Ax_l(dl_size, 0);
	SEQ_VECTOR <double> g_l(dl_size, 0);
	SEQ_VECTOR <double> Pg_l(dl_size, 0);
	SEQ_VECTOR <double> MPg_l(dl_size, 0);
	SEQ_VECTOR <double> z_l(dl_size, 0);
	SEQ_VECTOR <double> w_l(dl_size, 0);
	SEQ_VECTOR <double> Aw_l(dl_size, 0);
	SEQ_VECTOR <double> b_l(dl_size, 0);

	double gamma_l;
	double rho_l;
	double norm_l;
	double tol;
  double ztg = 0;
  double ztg_prew = 0;
  double wtAw;

	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );

  SparseMatrix W_l;
  W_l.type = 'G';
  W_l.rows = dl_size;
  W_l.cols = 0;


  SparseMatrix AW_l;
  AW_l.type = 'G';
  AW_l.rows = dl_size;
  AW_l.cols = 0;

	SEQ_VECTOR <double> Gamma_l  (dl_size, 0);
	SEQ_VECTOR <double> WtAW_l(dl_size, 0);

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d].norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ax_l[i] - b_l[i];
  }

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, g_l, Pg_l , 0);
	} else {
		Projector_l_compG    ( timeEvalProj, cluster, g_l, Pg_l , 0);
	}
	// *** Calculate the stop condition *******************************************
	tol = epsilon * parallel_norm_compressed(cluster, Pg_l);

	int precision = ceil(log(1 / epsilon) / log(10)) + 1;
	int iterationWidth = ceil(log(CG_max_iter) / log(10));
	std::string indent = "   ";

	auto spaces = [] (int count) {
		std::stringstream ss;
		for (int i = 0; i < count; i++) {
			ss << " ";
		}
		return ss.str();
	};

	ESINFO(CONVERGENCE)
		<< spaces(indent.size() + iterationWidth - 4) << "iter"
		<< spaces(indent.size() + precision - 3) << "|r|" << spaces(2)
		<< spaces(indent.size() + 4) << "r" << spaces(4)
		<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "e" << spaces(precision / 2)
		<< spaces(indent.size()) << "time[s]";

	// *** Start the CG iteration loop ********************************************
	for (int iter = -1; iter < CG_max_iter; iter++) {

		timing.totalTime.start();

    if (iter > -1) {
      W_l.dense_values.insert(W_l.dense_values.end(), w_l.begin(), w_l.end());
      W_l.nnz+=w_l.size();
      W_l.cols++;

      appA_time.start();
      apply_A_l_comp_dom_B(timeEvalAppa, cluster, w_l, Aw_l);
      appA_time.end();

      AW_l.dense_values.insert(AW_l.dense_values.end(), Aw_l.begin(), Aw_l.end());
      AW_l.nnz+=Aw_l.size();
      AW_l.cols++;

      wtAw = parallel_ddot_compressed(cluster, w_l, Aw_l);
      rho_l = -ztg/wtAw;

		  #pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
		  	x_l[i] = x_l[i] + rho_l * w_l[i];
        g_l[i] += Aw_l[i] * rho_l;
		  }
      ztg_prew = ztg;
    }
    switch (USE_PREC) {
    case config::solver::PRECONDITIONERalternative::LUMPED:
    case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
    case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
    case config::solver::PRECONDITIONERalternative::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, g_l, Pg_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, g_l, Pg_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      apply_prec_comp_dom_B(timeEvalPrec, cluster, Pg_l, MPg_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, MPg_l, z_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, MPg_l, z_l, 0 );
      }
      proj2_time.end();
      break;
    case config::solver::PRECONDITIONERalternative::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, g_l, z_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, g_l, z_l, 0 );
      }
      proj_time.end();
      break;
    default:
      ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
    }

    ztg = parallel_ddot_compressed(cluster, z_l, g_l);

    if (iter > -1) {
      gamma_l = ztg/ztg_prew;
		  #pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
		  	w_l[i] = z_l[i] +  w_l[i]*gamma_l;
		  }
    }
    else {
	    #pragma omp parallel for
for (size_t i = 0; i < w_l.size(); i++){
		  	w_l[i] = z_l[i];
		  }
    }

	  norm_time.start();
		norm_l = parallel_norm_compressed(cluster, Pg_l);
		norm_time.end();

		timing.totalTime.end();

		ESINFO(CONVERGENCE)
			<< indent << std::setw(iterationWidth) << iter + 1
			<< indent << std::fixed << std::setprecision(precision) <<  norm_l / tol * epsilon
			<< indent << std::scientific << std::setprecision(3) << norm_l
			<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
			<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();

		// *** Stop condition ******************************************************************
		if (norm_l < tol)
			break;

	} // end of CG iterations


	// *** save solution - in dual and amplitudes *********************************************


	#pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
		g_l[i] = -g_l[i];
	}


  dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = g_l;




	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case config::solver::PRECONDITIONERalternative::LUMPED:
	case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
	case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
	case config::solver::PRECONDITIONERalternative::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case config::solver::PRECONDITIONERalternative::NONE:
		timing.addEvent(proj_time);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

	timing.addEvent(appA_time );
	timing.addEvent(ddot_beta);
	timing.addEvent(ddot_alpha);

	// *** END - Preslocal out the timing for the iteration loop ***********************************

}
void IterSolverBase::Solve_full_ortho_CG_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                              C G      F U L L    O R T H O G O N A L                               # 
//##################################################################################################*/
//
	size_t dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);

	SEQ_VECTOR <double> Ax_l(dl_size, 0);
	SEQ_VECTOR <double> g_l(dl_size, 0);
	SEQ_VECTOR <double> Pg_l(dl_size, 0);
	SEQ_VECTOR <double> MPg_l(dl_size, 0);
	SEQ_VECTOR <double> z_l(dl_size, 0);
	SEQ_VECTOR <double> _z_l(dl_size, 0);
	SEQ_VECTOR <double> w_l(dl_size, 0);
	SEQ_VECTOR <double> Aw_l(dl_size, 0);
	SEQ_VECTOR <double> b_l(dl_size, 0);
	SEQ_VECTOR <double> v_tmp_l(dl_size, 0);

	SEQ_VECTOR <double> d_H(CG_max_iter, 0);
	SEQ_VECTOR <double> e_H(CG_max_iter, 0);



	double rho_l;
	double rho_l_prew = 1;
	double norm_l;
	double tol;
  double ztg;
  double wtAw;
  int cnt_iter=0;

	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );





  SparseMatrix W_l;
  W_l.type = 'G';
  W_l.rows = dl_size;
  W_l.cols = 0;


  SparseMatrix AW_l;
  AW_l.type = 'G';
  AW_l.rows = dl_size;
  AW_l.cols = 0;

	SEQ_VECTOR <double> Gamma_l  (CG_max_iter, 0);
	SEQ_VECTOR <double> _Gamma_l  (CG_max_iter, 0);
	SEQ_VECTOR <double> WtAW_l(CG_max_iter, 0);

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d].norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ax_l[i] - b_l[i];
  }

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, g_l, Pg_l , 0);
	} else {
		Projector_l_compG    ( timeEvalProj, cluster, g_l, Pg_l , 0);
	}
	// *** Calculate the stop condition *******************************************
	tol = epsilon * parallel_norm_compressed(cluster, Pg_l);

	int precision = ceil(log(1 / epsilon) / log(10)) + 1;
	int iterationWidth = ceil(log(CG_max_iter) / log(10));
	std::string indent = "   ";

	auto spaces = [] (int count) {
		std::stringstream ss;
		for (int i = 0; i < count; i++) {
			ss << " ";
		}
		return ss.str();
	};

	ESINFO(CONVERGENCE)
		<< spaces(indent.size() + iterationWidth - 4) << "iter"
		<< spaces(indent.size() + precision - 3) << "|r|" << spaces(2)
		<< spaces(indent.size() + 4) << "r" << spaces(4)
		<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "e" << spaces(precision / 2)
		<< spaces(indent.size()) << "time[s]";

	// *** Start the CG iteration loop ********************************************
	for (int iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.start();

    cnt_iter = iter - 1;

    if (iter > 0) {
      W_l.dense_values.insert(W_l.dense_values.end(), w_l.begin(), w_l.end());
      W_l.nnz+=w_l.size();
      W_l.cols++;

      appA_time.start();
      apply_A_l_comp_dom_B(timeEvalAppa, cluster, w_l, Aw_l);
      appA_time.end();

      AW_l.dense_values.insert(AW_l.dense_values.end(), Aw_l.begin(), Aw_l.end());
      AW_l.nnz+=Aw_l.size();
      AW_l.cols++;

      wtAw = parallel_ddot_compressed(cluster, w_l, Aw_l);
      WtAW_l[iter-1] = wtAw;

      ztg = parallel_ddot_compressed(cluster, z_l, g_l);


      rho_l = -ztg/wtAw;


      if (iter == 1)
      {
        d_H[iter-1] = -1.0/rho_l;
      }
      else
      {
        d_H[iter-1] = -(Gamma_l[iter-1]/rho_l_prew + 1.0/rho_l);
        e_H[iter-2] = sqrt(Gamma_l[iter-1])/rho_l_prew;
      }


      rho_l_prew = rho_l;


		  #pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
		  	x_l[i] = x_l[i] + rho_l * w_l[i];
        g_l[i] += Aw_l[i] * rho_l;
		  }
      //ztg_prew = ztg;
    }
    switch (USE_PREC) {
    case config::solver::PRECONDITIONERalternative::LUMPED:
    case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
    case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
    case config::solver::PRECONDITIONERalternative::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, g_l, Pg_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, g_l, Pg_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      apply_prec_comp_dom_B(timeEvalPrec, cluster, Pg_l, MPg_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, MPg_l, z_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, MPg_l, z_l, 0 );
      }
      proj2_time.end();
      break;
    case config::solver::PRECONDITIONERalternative::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, g_l, z_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, g_l, z_l, 0 );
      }
      Pg_l = z_l;
      proj_time.end();
      break;
    default:
      ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
    }

    if (iter > 0) {

      // filtering duplicit Lambda entries
      #pragma omp parallel for
for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++) {
        _z_l[i] = z_l[i] * cluster.my_lamdas_ddot_filter[i];
      }

      AW_l.DenseMatVec(_z_l,_Gamma_l,'T');

		  #pragma omp parallel for
for (eslocal i = 0; i < iter; i++) {
        _Gamma_l[i] /= -WtAW_l[i];
      }

	    MPI_Allreduce( &_Gamma_l[0], &Gamma_l[0], iter, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      W_l.DenseMatVec(Gamma_l,v_tmp_l);

		  #pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
		  	w_l[i] = z_l[i] +  v_tmp_l[i];
		  }

    }
    else {
	    #pragma omp parallel for
for (size_t i = 0; i < w_l.size(); i++){
		  	w_l[i] = z_l[i];
		  }
    }

	  norm_time.start();
		norm_l = parallel_norm_compressed(cluster, Pg_l);
		norm_time.end();

		timing.totalTime.end();

		ESINFO(CONVERGENCE)
			<< indent << std::setw(iterationWidth) << iter + 1
			<< indent << std::fixed << std::setprecision(precision) <<  norm_l / tol * epsilon
			<< indent << std::scientific << std::setprecision(3) << norm_l
			<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
			<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();

		// *** Stop condition ******************************************************************
		if (norm_l < tol)
			break;

	} // end of CG iterations


// EIGENVALUES AND EIGENVECTORS OF LANCZOS MATRIX
// Evaluation of cond(P*F*P) is limited by 1000 iter.
// Tridiagonal Lanczos' matrix is assembled at each node.
  bool cond_numb_FETI_operator=true;
  if (cnt_iter>0 && cnt_iter<1000 && cond_numb_FETI_operator && config::env::MPIrank==0){
    char JOBZ = 'N';
    double *Z = new double[cnt_iter];
    eslocal ldz = cnt_iter;
    LAPACKE_dstev(LAPACK_ROW_MAJOR, JOBZ, cnt_iter, &d_H[0], &e_H[0], Z, ldz);
    ESINFO(DETAILS) << "cond(P*F*P) = " << d_H[0]/d_H[cnt_iter-1]  ;
    delete [] Z;
  }


	// *** save solution - in dual and amplitudes *********************************************


	#pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
		g_l[i] = -g_l[i];
	}


  dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = g_l;




	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case config::solver::PRECONDITIONERalternative::LUMPED:
	case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
	case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
	case config::solver::PRECONDITIONERalternative::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case config::solver::PRECONDITIONERalternative::NONE:
		timing.addEvent(proj_time);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

	timing.addEvent(appA_time );
	timing.addEvent(ddot_beta);
	timing.addEvent(ddot_alpha);

	// *** END - Preslocal out the timing for the iteration loop ***********************************

}

void IterSolverBase::Solve_GMRES_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                                            G M R E S                                               # 
//##################################################################################################*/
//

	size_t dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);
	SEQ_VECTOR <double> Ax_l(dl_size, 0);
	SEQ_VECTOR <double> g_l(dl_size, 0);
	SEQ_VECTOR <double> Pg_l(dl_size, 0);
	SEQ_VECTOR <double> MPg_l(dl_size, 0);
	SEQ_VECTOR <double> Pw_l(dl_size, 0);
	SEQ_VECTOR <double> MPw_l(dl_size, 0);
	SEQ_VECTOR <double> z_l(dl_size, 0);
	SEQ_VECTOR <double> v_l(dl_size, 0);
	SEQ_VECTOR <double> w_l(dl_size, 0);
	SEQ_VECTOR <double> Aw_l(dl_size, 0);
        SEQ_VECTOR <double> b_l(dl_size, 0);
        SEQ_VECTOR <double> ones(dl_size, 1);

  int n_mat = CG_max_iter + 1;
	SEQ_VECTOR <double> b_H(n_mat, 0);
	SEQ_VECTOR <double> y_H(n_mat, 0);
	SEQ_VECTOR <double> g_H(n_mat, 0);

  double beta;
	double norm_l;
	double tol;
  double c_H,s_H;
  double norm_h;
  int cnt_iter=0;
  double tmp_double0;

  //cblas_
  double _alpha, _beta, _m, _n, _k, _lda, _ldb, _ldc;

	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );

#ifdef FLAG_VALIDATION
  SparseMatrix V_lt_V_l;
  V_lt_V_l.type = 'G';
  V_lt_V_l.rows = n_mat;
  V_lt_V_l.cols = n_mat;
  V_lt_V_l.dense_values.resize(n_mat*n_mat);
#endif

  SparseMatrix V_l;
  V_l.type = 'G';
  V_l.rows = dl_size;
  V_l.cols = 0;

  SEQ_VECTOR <double> Permut_l(n_mat*n_mat, 0);
  SEQ_VECTOR <double> Permut_tmp_l(n_mat*n_mat, 0);
  SEQ_VECTOR <double> P_tmp_P(n_mat*n_mat, 0);
  SEQ_VECTOR <double> H_l(n_mat*n_mat, 0);
  SEQ_VECTOR <double> H_l_modif(n_mat*n_mat, 0);
#ifdef FLAG_VALIDATION
  SEQ_VECTOR <double> tmp_H_l(n_mat*n_mat, 0);
#endif





    apply_A_l_comp_dom_B(timeEvalAppa, cluster, ones, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);



  for (int i = 0 ; i < n_mat; i++){
   Permut_l[n_mat*i + i] = 1;
   Permut_tmp_l[n_mat*i + i] = 1;
  }


	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d].norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ax_l[i] - b_l[i];
  }

  switch (USE_PREC) {
  case config::solver::PRECONDITIONERalternative::LUMPED:
  case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
  case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
  case config::solver::PRECONDITIONERalternative::MAGIC:
    proj1_time.start();
    if (USE_GGtINV == 1) {
      Projector_l_inv_compG( timeEvalProj, cluster, g_l, Pg_l, 0 );
    } else {
      Projector_l_compG		  ( timeEvalProj, cluster, g_l, Pg_l, 0 );
    }
    proj1_time.end();

    // Scale
    prec_time.start();
    apply_prec_comp_dom_B(timeEvalPrec, cluster, Pg_l, MPg_l);
    prec_time.end();
    // Re-Scale

    proj2_time.start();
    if (USE_GGtINV == 1) {
      Projector_l_inv_compG( timeEvalProj, cluster, MPg_l, z_l, 0 );
    } else {
      Projector_l_compG		  ( timeEvalProj, cluster, MPg_l, z_l, 0 );
    }
    proj2_time.end();
    break;
  case config::solver::PRECONDITIONERalternative::NONE:
    proj_time.start();
    if (USE_GGtINV == 1) {
      Projector_l_inv_compG( timeEvalProj, cluster, g_l, z_l, 0 );
    } else {
      Projector_l_compG		  ( timeEvalProj, cluster, g_l, z_l, 0 );
    }
    Pg_l = z_l;
    proj_time.end();
    break;
  default:
    ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
  }


	// *** Calculate the stop condition *******************************************
	//tol = epsilon * parallel_norm_compressed(cluster, Pg_l);

  norm_l = parallel_norm_compressed(cluster, z_l);
	tol = epsilon * norm_l;

	int precision = ceil(log(1 / epsilon) / log(10)) + 1;
	int iterationWidth = ceil(log(CG_max_iter) / log(10));
	std::string indent = "   ";

	auto spaces = [] (int count) {
		std::stringstream ss;
		for (int i = 0; i < count; i++) {
			ss << " ";
		}
		return ss.str();
	};

	ESINFO(CONVERGENCE)
		<< spaces(indent.size() + iterationWidth - 4) << "iter"
		<< spaces(indent.size() + precision - 3) << "|r|" << spaces(2)
		<< spaces(indent.size() + 4) << "r" << spaces(4)
		<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "e" << spaces(precision / 2)
		<< spaces(indent.size()) << "time[s]";



	//norm_l = parallel_norm_compressed(cluster, Pg_l);

  ESINFO(CONVERGENCE)
  	<< indent << std::setw(iterationWidth) << 1
  	<< indent << std::fixed << std::setprecision(precision) <<  1.0000000
  	<< indent << std::scientific << std::setprecision(3) << norm_l
  	<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
  	<< indent << std::fixed << std::setprecision(5) ;



  // initial gradient
  beta = sqrt(parallel_ddot_compressed(cluster, z_l, z_l));
  // InitialCondition for system H_{i+1,i} * y{i} = b_H
  b_H[0] = beta;
  // set-up first basis vector   (A * V_{i} = V_{i+1} * H_{i+1,i})
  #pragma omp parallel for
for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++) {
     v_l[i] = z_l[i]/beta;
  }



  tmp_double0 = parallel_ddot_compressed(cluster, v_l, v_l);


  V_l.dense_values.insert(V_l.dense_values.end(), v_l.begin(), v_l.end());
  V_l.nnz+=v_l.size();
  V_l.cols++;


//
  auto ij= [&]( eslocal ii, eslocal jj ) -> eslocal
   { return ii + n_mat*jj; };
  //
	// *** Start the CG iteration loop ********************************************
	for (int iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.start();

    appA_time.start();
    apply_A_l_comp_dom_B(timeEvalAppa, cluster, v_l, w_l);
    appA_time.end();

    switch (USE_PREC) {
    case config::solver::PRECONDITIONERalternative::LUMPED:
    case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
    case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
    case config::solver::PRECONDITIONERalternative::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, w_l, Pw_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, w_l, Pw_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      apply_prec_comp_dom_B(timeEvalPrec, cluster, Pw_l, MPw_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, MPw_l, z_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, MPw_l, z_l, 0 );
      }
      proj2_time.end();
      break;
    case config::solver::PRECONDITIONERalternative::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, w_l, z_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, w_l, z_l, 0 );
      }
      proj_time.end();
      break;
    default:
      ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
    }
//
//  Modified Gram-Schmidt
    for (int k = 0;k<iter+1;k++){
      H_l[ij(k,iter)] =parallel_ddot_compressed_double(cluster, &(V_l.dense_values[v_l.size()*k]), &(z_l[0]));

//
      #pragma omp parallel for
for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++) {
         z_l[i] -= V_l.dense_values[v_l.size()*k + i] * H_l[ij(k,iter)];

      }
    }
//
    H_l[ij(iter+1,iter)] = sqrt(parallel_ddot_compressed(cluster, z_l, z_l));
//
    #pragma omp parallel for
for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++) {
      v_l[i] = z_l[i]/H_l[ij(iter+1,iter)];
    }

    V_l.dense_values.insert(V_l.dense_values.end(), v_l.begin(), v_l.end());
    V_l.nnz+=v_l.size();
    V_l.cols++;

    // cblas set-up
     _alpha   = 1;
    _beta     = 0;
    _m        = iter+1;
    _n        = iter+1;
    _k        = iter+1;
    _lda      = n_mat;
    _ldb      = n_mat;
    _ldc      = n_mat;

    // next line isn't obligatory
    //
    w_l.insert(w_l.begin(),&(H_l[n_mat*iter]),&(H_l[n_mat*iter+iter+2]));
    //
		cblas_dgemv (CblasColMajor, CblasNoTrans, _m, _n,_alpha, &(Permut_l[0]), _lda, &(w_l[0]),
                                                  1,_beta, &(H_l_modif[n_mat*iter]), 1);
    //
    H_l_modif[n_mat*iter+iter+1] = H_l[n_mat*iter+iter+1];



    norm_h = sqrt(H_l_modif[ij(iter+1,iter)]*H_l_modif[ij(iter+1,iter)] +
                  H_l_modif[ij(iter,iter)]*H_l_modif[ij(iter,iter)]);
    c_H = H_l_modif[ij(iter  ,iter)] / norm_h;
    s_H = H_l_modif[ij(iter+1,iter)] / norm_h;


    if (iter>0){
      Permut_tmp_l[ij(iter-1,iter-1)] =  1;
      Permut_tmp_l[ij(iter  ,iter-1)] =  0;
      Permut_tmp_l[ij(iter-1,iter  )] =  0;
      Permut_tmp_l[ij(iter  ,iter  )] =  1;
    }

    Permut_tmp_l[ij(iter  ,iter  )] =  c_H;
    Permut_tmp_l[ij(iter+1,iter  )] = -s_H;
    Permut_tmp_l[ij(iter  ,iter+1)] =  s_H;
    Permut_tmp_l[ij(iter+1,iter+1)] =  c_H;

    H_l_modif[ij(iter  ,iter)] =  norm_h;
    H_l_modif[ij(iter+1,iter)] =  0;

    // cblas reset
    _m        = iter+2;
    _n        = iter+2;
    _k        = iter+2;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, _m, _n, _k, _alpha,
                      &(Permut_tmp_l[0]), _lda, &(Permut_l[0]), _ldb, _beta, &(P_tmp_P[0]), _ldc);

    // TODO: Do it better!
    for (int i = 0; i < n_mat*(iter+1)+iter+2;i++){
      Permut_l[i] = P_tmp_P[i];
    }

		cblas_dgemv (CblasColMajor, CblasNoTrans, _m, _n,_alpha, &(Permut_l[0]), _lda, &(b_H[0]),
                                                  1,_beta, &(g_H[0]), 1);


    norm_l = fabs(g_H[iter+1]);

		timing.totalTime.end();

		ESINFO(CONVERGENCE)
			<< indent << std::setw(iterationWidth) << iter + 2
			<< indent << std::fixed << std::setprecision(precision) <<  norm_l / tol * epsilon
			<< indent << std::scientific << std::setprecision(3) << norm_l
			<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
			<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();



#ifdef FLAG_VALIDATION
    _n = iter+1;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, _m, _n, _k, _alpha,
                      &(Permut_l[0]), _lda, &(H_l[0]), _ldb, _beta, &(tmp_H_l[0]), _ldc);
#endif

    cnt_iter = iter;
		if (norm_l < tol)
			break;
//
	}

  for (int i = cnt_iter-1;i>=0;i-- ){
    tmp_double0 = -g_H[i];
    for (int j = cnt_iter-1 ; j > i ; j-- ){
      tmp_double0-=H_l_modif[ij(i,j)] * y_H[j];
    }
    y_H[i] = tmp_double0 / H_l_modif[ij(i,i)];
  }

  _m = dl_size;
  _n = cnt_iter;
  _beta = 1;
  _lda = dl_size;
	cblas_dgemv (CblasColMajor, CblasNoTrans, _m, _n,_alpha, &(V_l.dense_values[0]), _lda, &(y_H[0]),
                                                  1,_beta, &(x_l[0]), 1);


  apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

  #pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
    w_l[i] = -(Ax_l[i] - b_l[i]);
  }



  dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = w_l;

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	}

#ifdef FLAG_VALIDATION
  for (int i = 0; i<cnt_iter+1;i++){
    for (int j = 0; j<cnt_iter+1;j++){
      w_l.insert(w_l.begin(),&(V_l.dense_values[dl_size*i]), &(V_l.dense_values[dl_size*(i+1)]));
      z_l.insert(  z_l.begin(),&(V_l.dense_values[dl_size*j]), &(V_l.dense_values[dl_size*(j+1)]));
      V_lt_V_l.dense_values[i*n_mat+j] = parallel_ddot_compressed(cluster, w_l, z_l);
    }
  }
#endif

  dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = w_l;


	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case config::solver::PRECONDITIONERalternative::LUMPED:
	case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
	case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
	case config::solver::PRECONDITIONERalternative::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case config::solver::PRECONDITIONERalternative::NONE:
		timing.addEvent(proj_time);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

	timing.addEvent(appA_time );
	timing.addEvent(ddot_beta);
	timing.addEvent(ddot_alpha);

	// *** END - Preslocal out the timing for the iteration loop ***********************************

} //  Solve_GMRES_singular_dom
//
void IterSolverBase::Solve_BICGSTAB_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                                         B I C G S T A B                                            # 
//##################################################################################################*/
//
//  FLAG_SOLUTION =  0;   solution found to the tolerance
//  FLAG_SOLUTION =  1;   no convergence for for given tolerance
//  FLAG_SOLUTION = -1;   breakdown: delta = 0
//  FLAG_SOLUTION = -2;   breakdown: omega = 0


	size_t dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);
	SEQ_VECTOR <double> Ay_l(dl_size, 0);
	SEQ_VECTOR <double> g_l(dl_size, 0);
	SEQ_VECTOR <double> Pg_l(dl_size, 0);
	SEQ_VECTOR <double> tmp1_l(dl_size, 0);
	SEQ_VECTOR <double> w_l(dl_size, 0);
	SEQ_VECTOR <double> t_l(dl_size, 0);
	SEQ_VECTOR <double> z_l(dl_size, 0);
	SEQ_VECTOR <double> v_l(dl_size, 0);
	SEQ_VECTOR <double> ztld_l(dl_size, 0);
	SEQ_VECTOR <double> b_l(dl_size, 0);
	SEQ_VECTOR <double> s_l(dl_size, 0);


	double norm_l; double tol;

	double rho      = 0;
	double omega    = 1;
	double delta    = 0;
	double delta_p  = 0;
	double gamma    = 0;
	int FLAG_SOLUTION = 0;


	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );

	if (USE_GGtINV == 1) {
		Projector_l_inv_compG( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector_l_compG	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ay_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d].norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ay_l[i] - b_l[i];
  }

  switch (USE_PREC) {
  case config::solver::PRECONDITIONERalternative::LUMPED:
  case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
  case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
  case config::solver::PRECONDITIONERalternative::MAGIC:
    proj1_time.start();
    if (USE_GGtINV == 1) {
      Projector_l_inv_compG( timeEvalProj, cluster, g_l, tmp1_l, 0 );
    } else {
      Projector_l_compG		  ( timeEvalProj, cluster, g_l, tmp1_l, 0 );
    }
    proj1_time.end();

    // Scale
    prec_time.start();
    apply_prec_comp_dom_B(timeEvalPrec, cluster, tmp1_l, g_l);
    prec_time.end();
    // Re-Scale

    proj2_time.start();
    if (USE_GGtINV == 1) {
      Projector_l_inv_compG( timeEvalProj, cluster, g_l, z_l, 0 );
    } else {
      Projector_l_compG		  ( timeEvalProj, cluster, g_l, z_l, 0 );
    }
    proj2_time.end();
    break;
  case config::solver::PRECONDITIONERalternative::NONE:
    proj_time.start();
    if (USE_GGtINV == 1) {
      Projector_l_inv_compG( timeEvalProj, cluster, g_l, z_l, 0 );
    } else {
      Projector_l_compG		  ( timeEvalProj, cluster, g_l, z_l, 0 );
    }
    Pg_l = z_l;
    proj_time.end();
    break;
  default:
    ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
  }


	// *** Calculate the stop condition *******************************************
	//tol = epsilon * parallel_norm_compressed(cluster, Pg_l);

  // norm_l = || b_bar ||
  norm_l = parallel_norm_compressed(cluster, z_l);

	tol = epsilon * norm_l;

	int precision = ceil(log(1 / epsilon) / log(10)) + 1;
	int iterationWidth = ceil(log(CG_max_iter) / log(10));
	std::string indent = "   ";

	auto spaces = [] (int count) {
		std::stringstream ss;
		for (int i = 0; i < count; i++) {
			ss << " ";
		}
		return ss.str();
	};

	ESINFO(CONVERGENCE)
		<< spaces(indent.size() + iterationWidth - 4) << "iter"
		<< spaces(indent.size() + precision - 3) << "|r|" << spaces(2)
		<< spaces(indent.size() + 4) << "r" << spaces(4)
		<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "e" << spaces(precision / 2)
		<< spaces(indent.size()) << "time[s]";



	//norm_l = parallel_norm_compressed(cluster, Pg_l);

  ESINFO(CONVERGENCE)
  	<< indent << std::setw(iterationWidth) << 1
  	<< indent << std::fixed << std::setprecision(precision) <<  1.0000000
  	<< indent << std::scientific << std::setprecision(3) << norm_l
  	<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
  	<< indent << std::fixed << std::setprecision(5) ;

  ztld_l = z_l;


  //
	// *** Start the CG iteration loop ********************************************
  int iter;
	for (iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.start();
    delta = parallel_ddot_compressed(cluster, z_l, ztld_l);

//    if (delta < 0.0001*epsilon){
    if (fabs(delta) == 1e-14){
      FLAG_SOLUTION=-1;
      break;
    }

    if (iter>0){
      gamma = -(delta/delta_p)*(rho/omega);
		  #pragma omp parallel for
for (size_t i = 0; i < w_l.size(); i++) {
        w_l[i] = z_l[i] + gamma * (w_l[i] - omega*v_l[i]);
      }
    }
    else
    {
      w_l = z_l;
    }

    appA_time.start();
    apply_A_l_comp_dom_B(timeEvalAppa, cluster, w_l, Ay_l);
    appA_time.end();

    switch (USE_PREC) {
    case config::solver::PRECONDITIONERalternative::LUMPED:
    case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
    case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
    case config::solver::PRECONDITIONERalternative::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, Ay_l, v_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, Ay_l, v_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      apply_prec_comp_dom_B(timeEvalPrec, cluster, v_l, tmp1_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, tmp1_l, v_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, tmp1_l, v_l, 0 );
      }
      proj2_time.end();
      break;
    case config::solver::PRECONDITIONERalternative::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, Ay_l, v_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, Ay_l, v_l, 0 );
      }
      proj_time.end();
      break;
    default:
      ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
    }

		timing.totalTime.end();


    rho = -(delta / parallel_ddot_compressed(cluster, ztld_l, v_l));


		#pragma omp parallel for
for (size_t i = 0; i < s_l.size(); i++) {
      s_l[i] = z_l[i] + rho * v_l[i];
    }

    norm_l = parallel_norm_compressed(cluster, s_l);
    
//  norm_l / tol * epsilon

    if (norm_l < tol){
      #pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
        x_l[i] += rho * w_l[i];
      } 
      //norm_l = parallel_norm_compressed(cluster, s_l);
      //error_ = norm_l / bnrm2;
      break;
    }
      
    appA_time.start();
    apply_A_l_comp_dom_B(timeEvalAppa, cluster, s_l, Ay_l);
    appA_time.end();
    
    switch (USE_PREC) {
    case config::solver::PRECONDITIONERalternative::LUMPED:
    case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
    case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
    case config::solver::PRECONDITIONERalternative::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, Ay_l, t_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, Ay_l, t_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      apply_prec_comp_dom_B(timeEvalPrec, cluster, t_l, tmp1_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, tmp1_l, t_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, tmp1_l, t_l, 0 );
      }
      proj2_time.end();
      break;
    case config::solver::PRECONDITIONERalternative::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_l_inv_compG( timeEvalProj, cluster, Ay_l, t_l, 0 );
      } else {
        Projector_l_compG		  ( timeEvalProj, cluster, Ay_l, t_l, 0 );
      }
      proj_time.end();
      break;
    default:
      ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
    }


    omega  = parallel_ddot_compressed(cluster, t_l, s_l)/parallel_ddot_compressed(cluster, t_l, t_l);


		#pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
      x_l[i] += rho * w_l[i] - omega * s_l[i];
      z_l[i] = -omega * t_l[i] + s_l[i];
    }
    
    norm_l = parallel_norm_compressed(cluster, z_l);



		ESINFO(CONVERGENCE)
			<< indent << std::setw(iterationWidth) << iter + 2
			<< indent << std::fixed << std::setprecision(precision) <<  norm_l / tol * epsilon
			<< indent << std::scientific << std::setprecision(3) << norm_l
			<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
			<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();


    
		if (norm_l< tol){
			break;
    }

    if ( fabs(omega) < 1e-14) {
      FLAG_SOLUTION=-2;
      break;
    }
    delta_p = delta;

//
	}


  if (iter == CG_max_iter){
      FLAG_SOLUTION=1;
  }


  ESINFO(CONVERGENCE)  << "FLAG_SOLUTION = " << FLAG_SOLUTION;
  


  dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = w_l;



  apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ay_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

  #pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
    w_l[i] = -(Ay_l[i] - b_l[i]);
  }



	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj , cluster, w_l, amplitudes, 2 );
	}


  dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = w_l;


	if (USE_GGtINV == 1) {
		Projector_l_inv_compG ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector_l_compG	  ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case config::solver::PRECONDITIONERalternative::LUMPED:
	case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
	case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
	case config::solver::PRECONDITIONERalternative::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case config::solver::PRECONDITIONERalternative::NONE:
		timing.addEvent(proj_time);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

	timing.addEvent(appA_time );
	timing.addEvent(ddot_beta);
	timing.addEvent(ddot_alpha);

	// *** END - Preslocal out the timing for the iteration loop ***********************************

} //  Solve_BICGSTAB_singular_dom



void IterSolverBase::Solve_PipeCG_singular_dom ( Cluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{

#ifdef USE_MPI_3
	if (mpi_rank == mpi_root)
		ESINFO(DETAILS) << "Note: PipeCG is using non-blocking AllReduce";
#endif
	size_t dl_size = cluster.my_lamdas_indices.size();

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
	double tol = 1;

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
	switch (USE_PREC) {
	case config::solver::PRECONDITIONERalternative::LUMPED:
	case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
	case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
	case config::solver::PRECONDITIONERalternative::MAGIC:

		#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
			tmp_l[i] = b_l[i] - Ax_l[i];
		}

		if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, r_l, 0 );
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, r_l, 0 );
		}

		tol = epsilon * parallel_norm_compressed(cluster, r_l);

		apply_prec_comp_dom_B(timeEvalPrec, cluster, r_l, tmp_l);
		if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, u_l, 0 );
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, u_l, 0 );
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, tmp_l); //apply_A_l_compB(timeEvalAppa, cluster, u_l, tmp_l);
		if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, w_l, 0 );
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, tmp_l, w_l, 0 );
		}

		break;
	case config::solver::PRECONDITIONERalternative::NONE:
		#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
			r_l[i] = b_l[i] - Ax_l[i];
		}

		if (USE_GGtINV == 1) {
			Projector_l_inv_compG( timeEvalProj, cluster, r_l, u_l, 0 );
		} else {
			Projector_l_compG    ( timeEvalProj, cluster, r_l, u_l, 0 );
		}
		tol = epsilon * parallel_norm_compressed(cluster, u_l);

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, w_l); 	//apply_A_l_compB(timeEvalAppa, cluster, u_l, w_l);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}


	int precision = ceil(log(1 / epsilon) / log(10)) + 1;
	int iterationWidth = ceil(log(CG_max_iter) / log(10));
	std::string indent = "   ";

	auto spaces = [] (int count) {
		std::stringstream ss;
		for (int i = 0; i < count; i++) {
			ss << " ";
		}
		return ss.str();
	};

	ESINFO(CONVERGENCE)
		<< spaces(indent.size() + iterationWidth - 4) << "iter"
		<< spaces(indent.size() + precision - 3) << "|r|" << spaces(2)
		<< spaces(indent.size() + 4) << "r" << spaces(4)
		<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "e" << spaces(precision / 2)
		<< spaces(indent.size()) << "time[s]";

	// *** Start the CG iteration loop ********************************************
	for (eslocal iter = 0; iter < CG_max_iter; iter++) {

		timing.totalTime.start();

		alpha_lp = alpha_l;
		gama_lp  = gama_l;

		//------------------------------------------
		ddot_time.start();
		MPI_Request mpi_req;

		SEQ_VECTOR <double> reduction_tmp (3,0);
		SEQ_VECTOR <double> send_buf      (3,0);

		switch (USE_PREC) {
		case config::solver::PRECONDITIONERalternative::LUMPED:
		case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
		case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
		case config::solver::PRECONDITIONERalternative::MAGIC:
			parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, r_l, &mpi_req, reduction_tmp, send_buf); // norm_l = parallel_norm_compressed(cluster, r_l);
			ddot_time.end();

			prec_time.start();
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, tmp_l);
			prec_time.end();

			proj_time.start();
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, m_l, 0 );
			} else {
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, m_l, 0 );
			}
			proj_time.end();

			appA_time.start();
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, tmp_l); //apply_A_l_compB(timeEvalAppa, cluster, m_l, tmp_l);
			appA_time.end();

			proj_time.start();
			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, tmp_l, n_l, 0 );
			} else {
				Projector_l_compG    ( timeEvalProj, cluster, tmp_l, n_l, 0 );
			}
			proj_time.end();

			break;
		case config::solver::PRECONDITIONERalternative::NONE:
			parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, u_l, &mpi_req, reduction_tmp, send_buf); // norm_l = parallel_norm_compressed(cluster, u_l);
			ddot_time.end();

			proj_time.start();

			if (USE_GGtINV == 1) {
				Projector_l_inv_compG( timeEvalProj, cluster, w_l, m_l, 0 );
			} else {
				Projector_l_compG    ( timeEvalProj, cluster, w_l, m_l, 0 );
			}

			proj_time.end();

			//------------------------------------------
			appA_time.start();
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, n_l); //apply_A_l_compB(timeEvalAppa, cluster, m_l, n_l);
			appA_time.end();
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
		}

#ifndef WIN32
#ifdef USE_MPI_3
		MPI_Wait(&mpi_req, &mpi_stat);
#endif
#endif

		norm_l  = sqrt(reduction_tmp[2]);
		if (norm_l < tol) {
			timing.totalTime.end();
			ESINFO(CONVERGENCE)
				<< indent << std::setw(iterationWidth) << iter + 1
				<< indent << std::fixed << std::setprecision(precision) <<  norm_l / tol * epsilon
				<< indent << std::scientific << std::setprecision(3) << norm_l
				<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
				<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();
			break;
		}

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

		#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
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

		norm_time.start();
		norm_time.end();

		 timing.totalTime.end();


		ESINFO(CONVERGENCE)
			<< indent << std::setw(iterationWidth) << iter + 1
			<< indent << std::fixed << std::setprecision(precision) <<  norm_l / tol * epsilon
			<< indent << std::scientific << std::setprecision(3) << norm_l
			<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
			<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();

	} // END of CG loop

	// *** save solution - in dual and amplitudes *********************************************

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l); //apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	#pragma omp parallel for
for(size_t i = 0; i < r_l.size(); i++)
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

	switch (USE_PREC) {
	case config::solver::PRECONDITIONERalternative::LUMPED:
	case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
	case config::solver::PRECONDITIONERalternative::DIRICHLET:
	case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
	case config::solver::PRECONDITIONERalternative::MAGIC:
		//timing.addEvent(proj1_time);
		timing.addEvent(proj_time);
		timing.addEvent(prec_time );
		//timing.addEvent(proj2_time);
		break;
	case config::solver::PRECONDITIONERalternative::NONE:
		timing.addEvent(proj_time);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

	timing.addEvent(appA_time);
	timing.addEvent(vec_time );

}


// *** Non-singular CG Solvers *******************************************
void IterSolverBase::Solve_RegCG_nonsingular  ( Cluster & cluster,
										    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
										    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel) {

	size_t dl_size = cluster.my_lamdas_indices.size();

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

	double beta_l  = 0;

	double alpha_l = 0;

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
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		// *** convert right hand side to dual
		cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[d]);

		y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
		cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
	}
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);



	eslocal iter = 0;
	timing.totalTime.reset();

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);

	#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {	// r = b - Ax;
		r_l[i] = b_l[i] - Ax_l[i];
		wp_l[i] = 0.0;
		yp_l[i] = 0.0;
	}

	tol = epsilon * parallel_norm_compressed(cluster, b_l);

	int precision = ceil(log(1 / epsilon) / log(10)) + 1;
	int iterationWidth = ceil(log(CG_max_iter) / log(10));
	std::string indent = "   ";

	auto spaces = [] (int count) {
		std::stringstream ss;
		for (int i = 0; i < count; i++) {
			ss << " ";
		}
		return ss.str();
	};

	ESINFO(CONVERGENCE)
		<< spaces(indent.size() + iterationWidth - 4) << "iter"
		<< spaces(indent.size() + 4) << "r" << spaces(4)
		<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "tol" << spaces(precision / 2)
		<< spaces(indent.size()) << "time[s]";

	for ( iter = 0; iter < 1000; iter++) {
		timing.totalTime.start();

		#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		switch (USE_PREC) {
		case config::solver::PRECONDITIONERalternative::LUMPED:
		case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
		case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
		case config::solver::PRECONDITIONERalternative::MAGIC:
			#pragma omp parallel for
for (size_t i = 0; i < w_l.size(); i++) {
				w_l[i] = r_l[i];
			}
			apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, y_l);
			break;
		case config::solver::PRECONDITIONERalternative::NONE:
			#pragma omp parallel for
for (size_t i = 0; i < w_l.size(); i++) {
				w_l[i] = r_l[i];
			}
			#pragma omp parallel for
for (size_t i = 0; i < w_l.size(); i++) {
				y_l[i] = w_l[i];
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
		}


		if (iter == 0) {									// if outputs.n_it==1;
			#pragma omp parallel for
for (size_t i = 0; i < y_l.size(); i++)
				p_l[i] = y_l[i];							// p = y;
		} else {
			ddot_beta.start();
			beta_l =          parallel_ddot_compressed(cluster, y_l, w_l);
			beta_l = beta_l / parallel_ddot_compressed(cluster, yp_l, wp_l);
			ddot_beta.end();

			#pragma omp parallel for
for (size_t i = 0; i < p_l.size(); i++)
				p_l[i] = y_l[i] + beta_l * p_l[i];			// p = y + beta * p;
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, p_l, Ap_l);



		alpha_l =           parallel_ddot_compressed(cluster, y_l, w_l);
		alpha_l = alpha_l / parallel_ddot_compressed(cluster, p_l, Ap_l);

		#pragma omp parallel for
for (size_t i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + alpha_l * p_l[i];
			r_l[i] = r_l[i] - alpha_l * Ap_l[i];
		}

		timing.totalTime.end();

		norm_l = parallel_norm_compressed(cluster, r_l);

		ESINFO(CONVERGENCE)
			<< indent << std::setw(iterationWidth) << iter + 1
			<< indent << std::scientific << std::setprecision(3) << norm_l
			<< indent << tol
			<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();


		if (norm_l < tol)
			break;

	} // end iter loop

	 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
	 timeGetSol.start();

	// reconstruction of u
	#pragma omp parallel for
for (size_t d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );

		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_l[ cluster.domains[d].lambda_map_sub_local[i]];

		cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');

		for(size_t i = 0; i < in_right_hand_side_primal[d].size(); i++)
			out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];

		cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);

	}

	 timeGetSol.endWithBarrier();
	 timing.addEvent(timeGetSol);


}

void IterSolverBase::Solve_PipeCG_nonsingular ( Cluster & cluster,
											SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
											SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel) {

		size_t dl_size = cluster.my_lamdas_indices.size();

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
//		for (size_t d = 1; d < cluster.domains.size(); d++) {
//			cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[0]);
//			cluster.domains[d].B1_comp.MatVec(cluster.x_prim_cluster1[0], cluster.compressed_tmp, 'N', 0, 0, 1.0);
//		}
//		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);

		//// *** convert right hand side to dual
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			// *** convert right hand side to dual
			cluster.domains[d].multKplusLocal(in_right_hand_side_primal[d],  cluster.x_prim_cluster1[d]);

			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, b_l);


		eslocal iter = 0;
		timing.totalTime.reset();

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);

		#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {	// r = b - Ax;
			r_l[i] = b_l[i] - Ax_l[i];
			wp_l[i] = 0.0;
			yp_l[i] = 0.0;
		}

		switch (USE_PREC) {
		case config::solver::PRECONDITIONERalternative::LUMPED:
		case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
		case config::solver::PRECONDITIONERalternative::DIRICHLET:
	  case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
		case config::solver::PRECONDITIONERalternative::MAGIC:
			apply_prec_comp_dom_B(timeEvalPrec, cluster, r_l, u_l);
			break;
		case config::solver::PRECONDITIONERalternative::NONE:
			#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
				u_l = r_l;
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, w_l);

		tol = epsilon * parallel_norm_compressed(cluster, b_l);

		int precision = ceil(log(1 / epsilon) / log(10)) + 1;
		int iterationWidth = ceil(log(CG_max_iter) / log(10));
		std::string indent = "   ";

		auto spaces = [] (int count) {
			std::stringstream ss;
			for (int i = 0; i < count; i++) {
				ss << " ";
			}
			return ss.str();
		};

		ESINFO(CONVERGENCE)
			<< spaces(indent.size() + iterationWidth - 4) << "iter"
			<< spaces(indent.size() + 4) << "r" << spaces(4)
			<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "tol" << spaces(precision / 2)
			<< spaces(indent.size()) << "time[s]";

		for ( iter = 0; iter < CG_max_iter; iter++) {

			timing.totalTime.start();

			alpha_lp = alpha_l;
			gama_lp  = gama_l;

			ddot_time.start();
			MPI_Request mpi_req;

			SEQ_VECTOR <double> reduction_tmp (2,0);
			SEQ_VECTOR <double> send_buf (2,0);
			parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, &mpi_req, reduction_tmp, send_buf);

			ddot_time.end();

			switch (USE_PREC) {
			case config::solver::PRECONDITIONERalternative::LUMPED:
			case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
			case config::solver::PRECONDITIONERalternative::DIRICHLET:
	    case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
			case config::solver::PRECONDITIONERalternative::MAGIC:
				prec_time.start();
				apply_prec_comp_dom_B(timeEvalPrec, cluster, w_l, m_l);
				prec_time.end();
				break;
			case config::solver::PRECONDITIONERalternative::NONE:
				#pragma omp parallel for
for (size_t i = 0; i < m_l.size(); i++) {
					m_l[i] = w_l[i];
				}
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
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

			#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
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

			ESINFO(CONVERGENCE)
				<< indent << std::setw(iterationWidth) << iter + 1
				<< indent << std::scientific << std::setprecision(3) << norm_l
				<< indent << tol
				<< indent << std::fixed << std::setprecision(5) << timing.totalTime.getLastStat();

			if (norm_l < tol)
				break;

		} // end iter loop

		// reconstruction of u

		 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
		 timeGetSol.start();

		#pragma omp parallel for
for (size_t d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );

			for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_l[ cluster.domains[d].lambda_map_sub_local[i]];

			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');

			for(size_t i = 0; i < in_right_hand_side_primal[d].size(); i++)
				out_primal_solution_parallel[d][i] = in_right_hand_side_primal[d][i] - cluster.x_prim_cluster1[d][i];

			cluster.domains[d].multKplusLocal(out_primal_solution_parallel[d]);

		}

		 timeGetSol.endWithBarrier();
		 timing.addEvent(timeGetSol);


}


// *** Coarse problem routines *******************************************
void IterSolverBase::CreateGGt( Cluster & cluster )

{
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
	eslocal count_cv = 0;
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

		ESINFO(PROGRESS2) << " Collecting matrices G : " << count_cv <<" of " << mpi_size;
	}

	//SparseMatrix Gtt;
	if (mpi_rank != mpi_root)
		G.Clear();
	//else {
	//Gtt = G;
	//G.Clear();
	//}
	// ****

	if (mpi_rank == mpi_root) {

		MKL_Set_Num_Threads(PAR_NUM_THREADS);
		// Create Gt and later GGt matrices and remove all elements under main diagonal of the GGt
		SparseMatrix Gt;

		double t1 = omp_get_wtime();
		G.MatTranspose(Gt);
		ESINFO(PROGRESS2) << "Gtranspose = " << omp_get_wtime() - t1;

		t1 = omp_get_wtime();
		SparseMatrix GGt_Mat;
		GGt_Mat.MatMat(G, 'N', Gt);
		ESINFO(PROGRESS2) << "G x Gt = " << omp_get_wtime() - t1;

		t1 = omp_get_wtime();
		Gt.Clear();
		G.Clear();
		ESINFO(PROGRESS2) << "G and Gt clear = " << omp_get_wtime() - t1;

		ESINFO(EXHAUSTIVE) << GGt_Mat.SpyText();

		t1 = omp_get_wtime();
		GGt_Mat.RemoveLower();
		ESINFO(PROGRESS2) << "GGt remove lower = " << omp_get_wtime() - t1;

		t1 = omp_get_wtime();
		// Create Sparse Direct solver for GGt
		GGt.msglvl = Info::report(LIBRARIES) ? 1 : 0;

		t1 = omp_get_wtime();
		GGt.ImportMatrix(GGt_Mat);
		ESINFO(PROGRESS2) << "ImportMatrix = " << omp_get_wtime() - t1;


		t1 = omp_get_wtime();
		GGt_Mat.Clear();


		t1 = omp_get_wtime();
		std::stringstream ss;
		ss << "Create GGt -> rank: " << config::env::MPIrank;
		GGt.Factorization(ss.str());
		ESINFO(PROGRESS2) << "Factorization = " << omp_get_wtime() - t1;


		t1 = omp_get_wtime();
		GGt.msglvl = 0;
		//TODO:
		MKL_Set_Num_Threads(1);
	}


	if (mpi_rank == mpi_root)
		GGtsize = GGt.cols;


	MPI_Bcast( & GGtsize, 1, esglobal_mpi, 0, MPI_COMM_WORLD);


#if TIME_MEAS >= 1
	double end = omp_get_wtime();
	ESINFO(PROGRESS2) <<"CG Loop - Create GGt  - collect all matrices   - Runtime = " << ec1 - sc1 << " s";
	ESINFO(PROGRESS2) <<"CG Loop - Create GGt  - GGt fact. processing   - Runtime = " << ep1 - sp1 << " s";
	ESINFO(PROGRESS2) <<"CG Loop - Create GGt  - total = proc + comm    - Runtime = " << end - start << " s";
#endif

}

void IterSolverBase::CreateGGt_inv_dist_d( Cluster & cluster )
{

	// temp variables
	vector < SparseMatrix > G_neighs   ( cluster.my_neighs.size() );
	vector < SparseMatrix > GGt_neighs ( cluster.my_neighs.size() );
	SparseMatrix Gt_l;
	SparseMatrix GGt_l;
	SparseMatrix GGt_Mat_tmp;
	SparseSolverCPU GGt_tmp;

    /* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs = Esutils::getEnv<int>("PAR_NUM_THREADS");
	GGt_tmp.iparm[2]  = num_procs;

	 TimeEvent SaRGlocal("Exchange local G1 matrices to neighs. "); SaRGlocal.start();
	ExchangeMatrices(cluster.G1, G_neighs, cluster.my_neighs);
	 SaRGlocal.end(); SaRGlocal.printStatMPI(); preproc_timing.addEvent(SaRGlocal);

	 TimeEvent Gt_l_trans("Local G1 matrix transpose to create Gt "); Gt_l_trans.start();
	if (cluster.USE_HFETI == 0) {
		cluster.G1.MatTranspose(Gt_l);
	}
	 Gt_l_trans.end(); Gt_l_trans.printStatMPI(); preproc_timing.addEvent(Gt_l_trans);

	 TimeEvent GxGtMatMat("Local G x Gt MatMat "); GxGtMatMat.start();
	if (cluster.USE_HFETI == 0) {
		GGt_l.MatMat(cluster.G1, 'N', Gt_l);
	} else {
		GGt_l.MatMatT(cluster.G1, cluster.G1);
	}
	 GxGtMatMat.end(); GxGtMatMat.printStatMPI(); preproc_timing.addEvent(GxGtMatMat);
	 //GxGtMatMat.PrintLastStatMPI_PerNode(0.0);

	for (size_t i = 0; i < GGt_l.CSR_J_col_indices.size(); i++) {
		GGt_l.CSR_J_col_indices[i] += mpi_rank * cluster.G1.rows;
	}
	GGt_l.cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;

	 TimeEvent GGTNeighTime("G1t_local x G1_neigh MatMat(N-times) "); GGTNeighTime.start();
	#pragma omp parallel for
for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {

		if (cluster.USE_HFETI == 0)
			GGt_neighs[neigh_i].MatMat(G_neighs[neigh_i], 'N', Gt_l);
		else
			GGt_neighs[neigh_i].MatMatT(G_neighs[neigh_i], cluster.G1);

		GGt_neighs[neigh_i].MatTranspose();

		eslocal inc = cluster.G1.rows * cluster.my_neighs[neigh_i];
		for (size_t i = 0; i < GGt_neighs[neigh_i].CSR_J_col_indices.size(); i++)
			GGt_neighs[neigh_i].CSR_J_col_indices[i] += inc;

		GGt_neighs[neigh_i].cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;
		G_neighs[neigh_i].Clear();
	}
	 GGTNeighTime.end(); GGTNeighTime.printStatMPI(); preproc_timing.addEvent(GGTNeighTime);
	 //GGTNeighTime.PrintLastStatMPI_PerNode(0.0);

	 TimeEvent GGtLocAsm("Assembling row of GGt per node - MatAddInPlace "); GGtLocAsm.start();
	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
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

		ESINFO(PROGRESS2) << "Collecting matrices G : " << count_cv_l <<" of " << mpi_size;
	}
	 collectGGt_time.end(); collectGGt_time.printStatMPI(); preproc_timing.addEvent(collectGGt_time);

	if (mpi_rank == 0)  {
		GGt_Mat_tmp.RemoveLower();
	}

	ESINFO(EXHAUSTIVE) << GGt_Mat_tmp.SpyText();

	MKL_Set_Num_Threads(PAR_NUM_THREADS);

	 TimeEvent GGt_bcast_time("Time to broadcast GGt from master all"); GGt_bcast_time.start();
	BcastMatrix(mpi_rank, mpi_root, mpi_root, GGt_Mat_tmp);
	 GGt_bcast_time.end(); GGt_bcast_time.printStatMPI(); preproc_timing.addEvent(GGt_bcast_time);

	// Create Sparse Direct solver for GGt
	if (mpi_rank == mpi_root) {
		GGt_tmp.msglvl = Info::report(LIBRARIES) ? 1 : 0;
	}

	 TimeEvent importGGt_time("Time to import GGt matrix into solver"); importGGt_time.start();
	GGt_tmp.ImportMatrix(GGt_Mat_tmp);
	 importGGt_time.end(); importGGt_time.printStatMPI(); preproc_timing.addEvent(importGGt_time);

	GGt_Mat_tmp.Clear();

	 TimeEvent GGtFactor_time("GGT Factorization time"); GGtFactor_time.start();
	 GGt_tmp.SetThreaded();
	 std::stringstream ss;
	 ss << "Create GGt_inv_dist-> rank: " << config::env::MPIrank;
	GGt_tmp.Factorization(ss.str());
	 GGtFactor_time.end();
	 //GGtFactor_time.printLastStatMPIPerNode();
	 GGtFactor_time.printStatMPI(); preproc_timing.addEvent(GGtFactor_time);

	 TimeEvent GGT_rhs_time("Time to create InitialCondition for get GGTINV"); GGT_rhs_time.start();
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

void IterSolverBase::CreateGGt_inv_dist( Cluster & cluster )
{

	// temp variables
	vector < SparseMatrix > G_neighs   ( cluster.my_neighs.size() );
	vector < SparseMatrix > GGt_neighs ( cluster.my_neighs.size() );
	SparseMatrix G1t_l;
	SparseMatrix GGt_l;
	SparseMatrix GGt_Mat_tmp;
	SparseSolverCPU GGt_tmp;

    /* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs = config::env::PAR_NUM_THREADS;
	GGt_tmp.iparm[2]  = num_procs;

	 TimeEvent SaRGlocal("Exchange local G1 matrices to neighs. "); SaRGlocal.start();
	if (cluster.SYMMETRIC_SYSTEM)  {
		ExchangeMatrices(cluster.G1, G_neighs, cluster.my_neighs);
	} else {
		ExchangeMatrices(cluster.G2, G_neighs, cluster.my_neighs);
	}
	 SaRGlocal.end(); SaRGlocal.printStatMPI(); preproc_timing.addEvent(SaRGlocal);

	 TimeEvent Gt_l_trans("Local G1 matrix transpose to create Gt "); Gt_l_trans.start();
	if (cluster.USE_HFETI == 0) {
		cluster.G1.MatTranspose(G1t_l);
	}
	 Gt_l_trans.end(); Gt_l_trans.printStatMPI(); preproc_timing.addEvent(Gt_l_trans);


	 if (cluster.SYMMETRIC_SYSTEM)  {
		  TimeEvent GxGtMatMat("Local G1 x G1t MatMat "); GxGtMatMat.start();
			if (cluster.USE_HFETI == 0) {
				GGt_l.MatMat(cluster.G1, 'N', G1t_l);
			} else {
				GGt_l.MatMatT(cluster.G1, cluster.G1);
			}
		  GxGtMatMat.end(); GxGtMatMat.printStatMPI(); preproc_timing.addEvent(GxGtMatMat);
	 } else {
		  TimeEvent GxGtMatMat("Local G2 x G1t MatMat "); GxGtMatMat.start();
			if (cluster.USE_HFETI == 0) {
				GGt_l.MatMat(cluster.G2, 'N', G1t_l);
			} else {
				GGt_l.MatMatT(cluster.G2, cluster.G1);
			}
                        GGt_l.MatTranspose(); 
		   GxGtMatMat.end(); GxGtMatMat.printStatMPI(); preproc_timing.addEvent(GxGtMatMat);
	 }
	 //GxGtMatMat.PrintLastStatMPI_PerNode(0.0);

	for (size_t i = 0; i < GGt_l.CSR_J_col_indices.size(); i++) {
		GGt_l.CSR_J_col_indices[i] += mpi_rank * cluster.G1.rows;
	}
	GGt_l.cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;

	 TimeEvent GGTNeighTime("G1t_local x G1_neigh MatMat(N-times) "); GGTNeighTime.start();
	#pragma omp parallel for
for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {

		GGt_neighs[neigh_i].MatMatT(G_neighs[neigh_i], cluster.G1);
		GGt_neighs[neigh_i].MatTranspose();

		//TODO: tady pocitam s tim, ze mam stejny pocet domen na cluster
		eslocal inc = cluster.G1.rows * cluster.my_neighs[neigh_i];
		for (size_t i = 0; i < GGt_neighs[neigh_i].CSR_J_col_indices.size(); i++)
			GGt_neighs[neigh_i].CSR_J_col_indices[i] += inc;

		//TODO: tady pocitam s tim, ze mam stejny pocet domen na cluster
		GGt_neighs[neigh_i].cols = cluster.NUMBER_OF_CLUSTERS * cluster.G1.rows;

		G_neighs[neigh_i].Clear();
	}
	 GGTNeighTime.end(); GGTNeighTime.printStatMPI(); preproc_timing.addEvent(GGTNeighTime);
	 //GGTNeighTime.PrintLastStatMPI_PerNode(0.0);

	 TimeEvent GGtLocAsm("Assembling row of GGt per node - MatAddInPlace "); GGtLocAsm.start();
	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		GGt_l.MatAddInPlace(GGt_neighs[neigh_i], 'N', 1.0);
		GGt_neighs[neigh_i].Clear();
	}
	 GGtLocAsm.end(); GGtLocAsm.printStatMPI(); preproc_timing.addEvent(GGtLocAsm);

	 TimeEvent collectGGt_time("Collect GGt pieces to master"); 	collectGGt_time.start();
	// Collecting pieces of GGt from all clusters to master (MPI rank 0) node - using binary tree reduction
	int count_cv_l = 0;
	for (eslocal li = 2; li <= 2*mpi_size; li = li * 2 ) {

		SparseMatrix recv_m_l;

		if (mpi_rank % li == 0) {
			if (li == 2) {
				GGt_Mat_tmp.MatAppend(GGt_l);
			}
			if ((mpi_rank + li/2) < mpi_size) {
				SendMatrix(mpi_rank, mpi_rank + li/2, GGt_l, mpi_rank,     recv_m_l);
				GGt_Mat_tmp.MatAppend(recv_m_l);
			} else {
				SendMatrix(mpi_rank, mpi_size + 1   , GGt_l, mpi_size + 1, recv_m_l);
			}
		} else {
			if ((mpi_rank + li/2) % li == 0) {
				if (li == 2) {
					SendMatrix(mpi_rank, mpi_rank       , GGt_l      , mpi_rank - li/2, recv_m_l);
				} else {
					SendMatrix(mpi_rank, mpi_rank       , GGt_Mat_tmp, mpi_rank - li/2, recv_m_l);
				}
			} else {
				SendMatrix(mpi_rank, mpi_rank+1, GGt_l, mpi_rank+1,recv_m_l);
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		GGt_l.Clear();

		count_cv_l += mpi_size/li;

		ESINFO(PROGRESS2) << "Collecting matrices G : " << count_cv_l <<" of " << mpi_size;
	}
	 collectGGt_time.end(); collectGGt_time.printStatMPI(); preproc_timing.addEvent(collectGGt_time);

	if (mpi_rank == 0 && cluster.SYMMETRIC_SYSTEM)  {
		ESINFO(EXHAUSTIVE) << "Creating symmetric Coarse problem (GGt) matrix";
		GGt_Mat_tmp.RemoveLower();
	} else {
		ESINFO(EXHAUSTIVE) << "Creating non-symmetric Coarse problem (GGt) matrix";
	}
	//Show GGt matrix structure in the solver LOG
	ESINFO(EXHAUSTIVE) << GGt_Mat_tmp.SpyText();

	// Entering data parallel region for single, in this case GGt matrix, we want MKL/Solver to run multi-threaded
	MKL_Set_Num_Threads(PAR_NUM_THREADS);

	//Broadcasting GGT matrix to all clusters/MPI ranks
	 TimeEvent GGt_bcast_time("Time to broadcast GGt from master all"); GGt_bcast_time.start();
	BcastMatrix(mpi_rank, mpi_root, mpi_root, GGt_Mat_tmp);
	 GGt_bcast_time.end(); GGt_bcast_time.printStatMPI(); preproc_timing.addEvent(GGt_bcast_time);


	// *** Calculating inverse GGt matrix in distributed fashion ***********************************************************
	// Create Sparse Direct solver for GGt
	if (mpi_rank == mpi_root) {
		GGt_tmp.msglvl = Info::report(LIBRARIES) ? 1 : 0;
	}

	 TimeEvent importGGt_time("Time to import GGt matrix into solver"); importGGt_time.start();
	//TODO: Toto se musi nejak korektne vyresit pomoci fyziky
//	if (cluster.SYMMETRIC_SYSTEM)  {
//		GGt_tmp.mtype = 2; // non-symmetric; // GGt_tmp.mtype = 2;  // Real symmetric positive definite matrix - this is default setting of the PARDISO solver
//		GGt_Mat_tmp.mtype = SparseMatrix::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;
//	} else {
//		GGt_tmp.mtype = 11; // non-symmetric
//		GGt_Mat_tmp.mtype = SparseMatrix::MatrixType::REAL_UNSYMMETRIC;
//	}
	GGt_Mat_tmp.mtype = cluster.mtype;

	GGt_tmp.ImportMatrix_wo_Copy (GGt_Mat_tmp);
	 importGGt_time.end(); importGGt_time.printStatMPI(); preproc_timing.addEvent(importGGt_time);

	//GGt_Mat_tmp.Clear();

	 TimeEvent GGtFactor_time("GGT Factorization time"); GGtFactor_time.start();
	 GGt_tmp.SetThreaded();
	 std::stringstream ss;
	 ss << "Create GGt_inv_dist-> rank: " << config::env::MPIrank;
	GGt_tmp.Factorization(ss.str());
	 GGtFactor_time.end();
	 //GGtFactor_time.printLastStatMPIPerNode();
	 GGtFactor_time.printStatMPI(); preproc_timing.addEvent(GGtFactor_time);

	 TimeEvent GGT_rhs_time("Time to create InitialCondition for get GGTINV"); GGT_rhs_time.start();
	 //TODO: tady pocitam s tim, ze mam stejny pocet domen na cluster
	SEQ_VECTOR <double> rhs   (cluster.G1.rows * GGt_tmp.rows, 0.0);
	cluster.GGtinvV.resize    (cluster.G1.rows * GGt_tmp.rows, 0.0);

	//TODO: tady pocitam s tim, ze mam stejny pocet domen na cluster
	for (eslocal i = 0; i < cluster.G1.rows; i++) {
		eslocal index = (GGt_tmp.rows * i) + (cluster.G1.rows * mpi_rank) + i;
		rhs[index] = 1;
	}
	GGT_rhs_time.end(); GGT_rhs_time.printStatMPI(); preproc_timing.addEvent(GGT_rhs_time);

	 TimeEvent GGt_solve_time("Running solve to get stripe(s) of GGtINV"); GGt_solve_time.start();

	GGt_tmp.Solve(rhs, cluster.GGtinvV, cluster.G1.rows);

	cluster.GGtinvM.dense_values = cluster.GGtinvV;
	cluster.GGtinvM.cols 		 = cluster.G1.rows;
	cluster.GGtinvM.rows	 	 = GGt_tmp.rows;
    cluster.GGtinvM.type 		 = 'G';

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
void IterSolverBase::Projector_l_compG (TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
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
			#pragma omp parallel for
for (size_t i = 0; i < x_in.size(); i++)
				y_out[i] = x_in[i] - y_out[i];
		}

	}

	time_eval.totalTime.end();

}

void IterSolverBase::Projector_l_inv_compG (TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
{

	 time_eval.totalTime.start();

	eslocal d_local_size = cluster.G1_comp.rows;

	SEQ_VECTOR<double> d_local( d_local_size );
	SEQ_VECTOR<double> d_mpi  ( GGtsize );

	 time_eval.timeEvents[0].start();
	if ( output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 1 || output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 3) {
		d_local = x_in;
	} else {
		if (cluster.SYMMETRIC_SYSTEM) {
			cluster.G1_comp.MatVec(x_in, d_local, 'N');
		} else {
			cluster.G2_comp.MatVec(x_in, d_local, 'N');
		}
	}
	 time_eval.timeEvents[0].end();

	//TODO: Pocitam s tim, ze kazdy cluster ma stejny ocet domen
	 time_eval.timeEvents[1].start();
	MPI_Allgather(&d_local[0], d_local_size, MPI_DOUBLE,
		&d_mpi[0], d_local_size, MPI_DOUBLE,
		MPI_COMM_WORLD);
	 time_eval.timeEvents[1].end();

	 time_eval.timeEvents[2].start();
	cluster.GGtinvM.DenseMatVec(d_mpi, d_local, 'T');
	 time_eval.timeEvents[2].end();

	 time_eval.timeEvents[3].start();
	//MPI_Scatter( &d_mpi[0],      d_local_size, MPI_DOUBLE,
	//	&d_local[0], d_local_size, MPI_DOUBLE,
	//	mpi_root, MPI_COMM_WORLD);
	 time_eval.timeEvents[3].end();

	if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 2 || output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 3) {
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
			#pragma omp parallel for
for (size_t i = 0; i < y_out.size(); i++)
				y_out[i] = x_in[i] - y_out[i];
		}

	}

	time_eval.totalTime.end();
}


void IterSolverBase::Projector_l_inv_compG_d (TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
{

	time_eval.totalTime.start();

	eslocal d_local_size = cluster.G1_comp.rows;

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
			#pragma omp parallel for
for (size_t i = 0; i < y_out.size(); i++)
				y_out[i] = x_in[i] - y_out[i];
		}

	}

	time_eval.totalTime.end();

}



void IterSolverBase::apply_prec_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ) {

	time_eval.totalTime.start();

	time_eval.timeEvents[0].start();

	#pragma omp parallel for
for (size_t d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling

		switch (USE_PREC) {
		case config::solver::PRECONDITIONERalternative::LUMPED:
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
			cluster.domains[d].K.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			cluster.domains[d]._RegMat.MatVecCOO(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N', -1.0);
			break;
		case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster2[d], 'T');
			break;
		case config::solver::PRECONDITIONERalternative::DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
			//cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			cluster.domains[d].Prec.DenseMatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			break;
		case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
			cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			break;
		case config::solver::PRECONDITIONERalternative::MAGIC:
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
			cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			break;
		case config::solver::PRECONDITIONERalternative::NONE:
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
		}

	}

	std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
	SEQ_VECTOR < double > y_out_tmp;
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );


		switch (USE_PREC) {
		case config::solver::PRECONDITIONERalternative::LUMPED:
		case config::solver::PRECONDITIONERalternative::WEIGHT_FUNCTION:
		case config::solver::PRECONDITIONERalternative::MAGIC:
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
    //TODO  check if MatVec is correct (DenseMatVec!!!) 
		case config::solver::PRECONDITIONERalternative::DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
		case config::solver::PRECONDITIONERalternative::SUPER_DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
		case config::solver::PRECONDITIONERalternative::NONE:
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
		}


		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
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

namespace espreso {


// **********************************************************************
// *** Communication layer **********************************************
void   SendMatrix2  ( eslocal rank, eslocal source_rank, SparseMatrix & A_in, eslocal dest_rank, SparseMatrix & B_out) {

	eslocal param_tag = 1;
	eslocal I_row_tag = 2;
	eslocal J_col_tag = 3;
	eslocal V_val_tag = 4;

	MPI_Status status;

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

void   SendMatrix  ( eslocal rank, eslocal source_rank, SparseMatrix & A_in, eslocal dest_rank, SparseMatrix & B_out) {

	eslocal tag = 1;

	if (rank == source_rank) {

		SEQ_VECTOR < MPI_Request > request ( 4 );
		eslocal send_par_buf[4];

		send_par_buf[0] = A_in.cols;
		send_par_buf[1] = A_in.rows;
		send_par_buf[2] = A_in.nnz;
		send_par_buf[3] = A_in.type;

		MPI_Isend(send_par_buf, 		   				  4, 	esglobal_mpi, 	dest_rank, tag, MPI_COMM_WORLD, &request[0] );
		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, 	esglobal_mpi, 	dest_rank, tag, MPI_COMM_WORLD, &request[1] );
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      	esglobal_mpi, 	dest_rank, tag, MPI_COMM_WORLD, &request[2] );
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   		MPI_DOUBLE, 	dest_rank, tag, MPI_COMM_WORLD, &request[3] );

		MPI_Waitall( 4 , &request[0], MPI_STATUSES_IGNORE);

	}

	if (rank == dest_rank) {

		MPI_Status status;
		SEQ_VECTOR < MPI_Request > request ( 3 );
		eslocal recv_par_buf[4];

		MPI_Recv(recv_par_buf, 4, esglobal_mpi, source_rank, tag, MPI_COMM_WORLD, & status);
		B_out.cols = recv_par_buf[0];
		B_out.rows = recv_par_buf[1];
		B_out.nnz  = recv_par_buf[2];
		B_out.type = recv_par_buf[3];

		B_out.CSR_I_row_indices.resize(B_out.rows + 1);
		B_out.CSR_J_col_indices.resize(B_out.nnz);
		B_out.CSR_V_values.     resize(B_out.nnz);

		MPI_Irecv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, esglobal_mpi,    source_rank, tag, MPI_COMM_WORLD, &request[0] );
		MPI_Irecv(&B_out.CSR_J_col_indices[0], B_out.nnz,      esglobal_mpi,    source_rank, tag, MPI_COMM_WORLD, &request[1] );
		MPI_Irecv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE,      source_rank, tag, MPI_COMM_WORLD, &request[2] );

		MPI_Waitall( 3 , &request[0], MPI_STATUSES_IGNORE);

	}

}


void   SendMatrix ( SparseMatrix & A_in, eslocal dest_rank ) {

	int rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);


	eslocal param_tag = 1;
	eslocal I_row_tag = 2;
	eslocal J_col_tag = 3;
	eslocal V_val_tag = 4;

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

void ExchangeMatrices (SparseMatrix & A_in, SEQ_VECTOR <SparseMatrix> & B_out, SEQ_VECTOR <eslocal> neighbor_ranks ) {

	eslocal tag = 1;

	SEQ_VECTOR < MPI_Request > request ( 7 * neighbor_ranks.size() );

	SEQ_VECTOR < SEQ_VECTOR < eslocal > > send_par_buf ( neighbor_ranks.size() );
	SEQ_VECTOR < SEQ_VECTOR < eslocal > > recv_par_buf ( neighbor_ranks.size() );

	//Send Matrix properties
	for (size_t neigh_i = 0; neigh_i < neighbor_ranks.size(); neigh_i++ )
	{
		send_par_buf[neigh_i].resize(4);
		eslocal dest_rank = neighbor_ranks[neigh_i];

		send_par_buf[neigh_i][0] = A_in.cols;
		send_par_buf[neigh_i][1] = A_in.rows;
		send_par_buf[neigh_i][2] = A_in.nnz;
		send_par_buf[neigh_i][3] = (eslocal)A_in.type;

		MPI_Isend(&send_par_buf[neigh_i][0],  4, 				esglobal_mpi, 	dest_rank, tag, MPI_COMM_WORLD, &request[7 * neigh_i + 0] );

	}


	for (size_t neigh_i = 0; neigh_i < neighbor_ranks.size(); neigh_i++ )
	{
		recv_par_buf[neigh_i].resize(4);
		eslocal source_rank = neighbor_ranks[neigh_i];

		MPI_Recv(&recv_par_buf[neigh_i][0], 4, esglobal_mpi, source_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		B_out[neigh_i].cols = recv_par_buf[neigh_i][0];
		B_out[neigh_i].rows = recv_par_buf[neigh_i][1];
		B_out[neigh_i].nnz  = recv_par_buf[neigh_i][2];
		B_out[neigh_i].type = (char)recv_par_buf[neigh_i][3];

		B_out[neigh_i].CSR_I_row_indices.resize(B_out[neigh_i].rows + 1);
		B_out[neigh_i].CSR_J_col_indices.resize(B_out[neigh_i].nnz);
		B_out[neigh_i].CSR_V_values.     resize(B_out[neigh_i].nnz);
	}

	// Send Data
	for (size_t neigh_i = 0; neigh_i < neighbor_ranks.size(); neigh_i++ )
	{
		eslocal dest_rank = neighbor_ranks[neigh_i];
		eslocal tag = 1;

		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, 	esglobal_mpi, 	dest_rank, tag, MPI_COMM_WORLD, &request[7 * neigh_i + 1] );
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      	esglobal_mpi, 	dest_rank, tag, MPI_COMM_WORLD, &request[7 * neigh_i + 2] );
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   		MPI_DOUBLE, 	dest_rank, tag, MPI_COMM_WORLD, &request[7 * neigh_i + 3] );

	}


	for (size_t neigh_i = 0; neigh_i < neighbor_ranks.size(); neigh_i++ )
	{
		eslocal source_rank = neighbor_ranks[neigh_i];
		eslocal tag = 1;

		MPI_Irecv(&B_out[neigh_i].CSR_I_row_indices[0], B_out[neigh_i].rows + 1, esglobal_mpi,    source_rank, tag, MPI_COMM_WORLD, &request[7 * neigh_i + 4] );
		MPI_Irecv(&B_out[neigh_i].CSR_J_col_indices[0], B_out[neigh_i].nnz,      esglobal_mpi,    source_rank, tag, MPI_COMM_WORLD, &request[7 * neigh_i + 5] );
		MPI_Irecv(&B_out[neigh_i].CSR_V_values[0],      B_out[neigh_i].nnz,      MPI_DOUBLE,      source_rank, tag, MPI_COMM_WORLD, &request[7 * neigh_i + 6] );
	}

	MPI_Waitall(7 * neighbor_ranks.size(), &request[0], MPI_STATUSES_IGNORE);

}

void   BcastMatrix ( eslocal rank, eslocal mpi_root, eslocal source_rank, SparseMatrix & A) {

	eslocal send_par_buf[4];

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

void   All_Reduce_lambdas_compB2( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out )
{

	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			cluster.my_comm_lambdas[i][j] = x_in[cluster.my_comm_lambdas_indices_comp[i][j]];
		}
	}


	MPI_Request * mpi_req  = new MPI_Request [cluster.my_neighs.size()];
	MPI_Status  * mpi_stat = new MPI_Status  [cluster.my_neighs.size()];

	cluster.iter_cnt_comm++;
	eslocal tag = cluster.iter_cnt_comm;

	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		MPI_Sendrecv(
			&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,
			&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,
			MPI_COMM_WORLD, &mpi_stat[neigh_i] );
	}

	//for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
	//	size_t b_size = cluster.my_comm_lambdas[neigh_i].size();
	//	MPI_Isend(&b_size,                              1                                      , esglobal_mpi   , cluster.my_neighs[neigh_i], tag + 100, MPI_COMM_WORLD, &mpi_req[neigh_i] );
	//	MPI_Isend(&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,       MPI_COMM_WORLD, &mpi_req[neigh_i] );
	//
	//}

	//for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
	//	size_t r_size = 0;
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
	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			y_out[cluster.my_comm_lambdas_indices_comp[i][j]] += cluster.my_recv_lambdas[i][j];
		}
	}



}


void   All_Reduce_lambdas_compB( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out )
{

	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			cluster.my_comm_lambdas[i][j] = x_in[cluster.my_comm_lambdas_indices_comp[i][j]];
		}
	}

	SEQ_VECTOR < MPI_Request > request ( 2 * cluster.my_neighs.size() );

	cluster.iter_cnt_comm++;
	eslocal tag = 1;

	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		MPI_Isend(
			&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag, MPI_COMM_WORLD, &request[ 0                        + neigh_i] );
	}

	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		MPI_Irecv(
			&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag, MPI_COMM_WORLD, &request[ cluster.my_neighs.size() + neigh_i] );
	}

	MPI_Waitall( 2 * cluster.my_neighs.size(), &request[0], MPI_STATUSES_IGNORE);

	y_out = x_in; // TODO: POZOR pozor
	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			y_out[cluster.my_comm_lambdas_indices_comp[i][j]] += cluster.my_recv_lambdas[i][j];
		}
	}

}

void   compress_lambda_vector  ( Cluster & cluster, SEQ_VECTOR <double> & decompressed_vec_lambda)
{
	//compress vector for CG in main loop
	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[cluster.my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(cluster.my_lamdas_indices.size());
}

void   decompress_lambda_vector( Cluster & cluster, SEQ_VECTOR <double> & compressed_vec_lambda)
{
	SEQ_VECTOR <double> decompressed_vec_lambda (cluster.domains[0].B1.rows,0);

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[cluster.my_lamdas_indices[i]] = compressed_vec_lambda[i];

	compressed_vec_lambda = decompressed_vec_lambda;
}

double parallel_norm_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector )
{

	double wl = 0; double wg = 0;

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)
		wl = wl + (input_vector[i] * input_vector[i] * cluster.my_lamdas_ddot_filter[i]);

	MPI_Allreduce( &wl, &wg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double norm_l = sqrt(wg);

	return norm_l;
}


double parallel_ddot_compressed_double( Cluster & cluster, double * input_vector1, double * input_vector2 )
{
	double a1 = 0; double a1g = 0;

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
		a1 = a1 + (input_vector1[i] * input_vector2[i] * cluster.my_lamdas_ddot_filter[i]);
	}

	MPI_Allreduce( &a1, &a1g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return a1g;
}


double parallel_ddot_compressed( Cluster & cluster, SEQ_VECTOR<double> & input_vector1, SEQ_VECTOR<double> & input_vector2 )
{
	double a1 = 0; double a1g = 0;

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
		a1 = a1 + (input_vector1[i] * input_vector2[i] * cluster.my_lamdas_ddot_filter[i]);
	}

	MPI_Allreduce( &a1, &a1g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return a1g;
}

void   parallel_ddot_compressed_non_blocking( Cluster & cluster,
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b,
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b,
	SEQ_VECTOR<double> & input_norm_vec,

	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf)
{

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
		send_buf[0] = send_buf[0] + (input_vector_1a[i] * input_vector_1b[i] * cluster.my_lamdas_ddot_filter[i]); // ddot 1
		send_buf[1] = send_buf[1] + (input_vector_2a[i] * input_vector_2b[i] * cluster.my_lamdas_ddot_filter[i]); // ddot 2
		send_buf[2] = send_buf[2] + (input_norm_vec[i]  * input_norm_vec[i]  * cluster.my_lamdas_ddot_filter[i]); // norm
	}


#ifdef WIN32
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce( &send_buf[0], &output[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
#ifdef USE_MPI_3
	MPI_Iallreduce( &send_buf[0], &output[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, mpi_req);
#else
	MPI_Allreduce( &send_buf[0], &output[0], 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif

}


void   parallel_ddot_compressed_non_blocking( Cluster & cluster,
	SEQ_VECTOR<double> & input_vector_1a, SEQ_VECTOR<double> & input_vector_1b,
	SEQ_VECTOR<double> & input_vector_2a, SEQ_VECTOR<double> & input_vector_2b,

	MPI_Request * mpi_req,
	SEQ_VECTOR <double> & output,
	SEQ_VECTOR <double> & send_buf)
{

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
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

}

// *** END - Communication layer ****************************************
// **********************************************************************





