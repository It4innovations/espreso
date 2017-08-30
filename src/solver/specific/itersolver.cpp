//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

#include "itersolver.h"

#include "../../basis/utilities/utils.h"

#ifdef READEX_LEVEL_1
#include "readex.h"
#include "readex_regions.h"
#endif

using namespace espreso;

IterSolverBase::IterSolverBase(const ESPRESOSolver &configuration):
	configuration(configuration),
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

	MPI_Comm_rank (environment->MPICommunicator, &mpi_rank);	/* get current process id */
	MPI_Comm_size (environment->MPICommunicator, &mpi_size);	/* get number of processes */
	mpi_root = 0;

}




void IterSolverBase::Preprocessing ( SuperCluster & cluster )
{
	// Coarse problem - Make GGt

	 preproc_timing.totalTime.start();
	 TimeEvent createGGT_time("Time to create GGt");createGGT_time.start();

	if (USE_GGtINV == 1) {
		CreateGGt_Inv( cluster );
	} else {
		ESINFO(GLOBAL_ERROR) << "Only Inverse of GGT is supported for Projector";
		// CreateGGt    ( cluster );
	}

	 createGGT_time.end();
	 createGGT_time.printStatMPI();
	 preproc_timing.addEvent(createGGT_time);
	 preproc_timing.totalTime.end();


}



void IterSolverBase::Solve ( SuperCluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & out_primal_solution_parallel,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & out_dual_solution_parallel)
{

#ifdef READEX_LEVEL_1
	READEX_REGION_START( REG_Solve_RegCG, "Solve_RegCG", SCOREP_USER_REGION_TYPE_COMMON );
#endif

	switch (configuration.solver) {
	case ESPRESO_ITERATIVE_SOLVER::PCG:
		Solve_RegCG ( cluster, in_right_hand_side_primal );
		break;
	case ESPRESO_ITERATIVE_SOLVER::pipePCG:
		Solve_PipeCG_singular_dom( cluster, in_right_hand_side_primal );
		break;
	case ESPRESO_ITERATIVE_SOLVER::orthogonalPCG:
		Solve_full_ortho_CG_singular_dom (cluster, in_right_hand_side_primal );
		break;
	case ESPRESO_ITERATIVE_SOLVER::GMRES:
		Solve_GMRES_singular_dom (cluster, in_right_hand_side_primal );
		break;
	case ESPRESO_ITERATIVE_SOLVER::BICGSTAB:
		Solve_BICGSTAB_singular_dom(cluster, in_right_hand_side_primal );
		break;
	case ESPRESO_ITERATIVE_SOLVER::QPCE:
		Solve_QPCE_singular_dom(cluster, in_right_hand_side_primal );
		break;
	case ESPRESO_ITERATIVE_SOLVER::orthogonalPCG_CP:
		Solve_full_ortho_CG_singular_dom_geneo(cluster, in_right_hand_side_primal);
		break;
	case ESPRESO_ITERATIVE_SOLVER::PCG_CP:
		ESINFO(GLOBAL_ERROR) << "Regular CG with conjugate projector not implemented yet";
		break;

	default:
		ESINFO(GLOBAL_ERROR) << "Unknown CG solver";
	}

	 postproc_timing.totalTime.start();

	 TimeEvent timeGetSol(string("Solver - Get Primal Solution"));
	 timeGetSol.start();
	GetSolution_Primal_singular_parallel( cluster, in_right_hand_side_primal, out_primal_solution_parallel, out_dual_solution_parallel );
	 timeGetSol.endWithBarrier();
	 postproc_timing.addEvent(timeGetSol);

	 postproc_timing.totalTime.endWithBarrier();

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Solve_RegCG);
#endif

}


void IterSolverBase::GetResiduum_Dual_singular_parallel    ( SuperCluster & cluster, SEQ_VECTOR <double> & dual_residuum_out ) {

	dual_residuum_out = dual_residuum_compressed_parallel;
	cluster.decompress_lambda_vector( dual_residuum_out );

}

void IterSolverBase::GetSolution_Dual_singular_parallel    ( SuperCluster & cluster, SEQ_VECTOR <double> & dual_solution_out, SEQ_VECTOR<double> & amplitudes_out ) {

	cluster.decompress_lambda_vector( dual_solution_out );
	dual_solution_out = dual_soultion_compressed_parallel;

	amplitudes_out	  = amplitudes;

}

void IterSolverBase::GetSolution_Primal_singular_parallel  ( SuperCluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
		SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out,
		SEQ_VECTOR < SEQ_VECTOR <double> > & dual_solution_out) {

	MakeSolution_Primal_singular_parallel(cluster, in_right_hand_side_primal, primal_solution_out );

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
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d]->B1_comp_dom.rows, 0.0 );
			for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = dual_soultion_compressed_parallel[ cluster.domains[d]->lambda_map_sub_local[i]]; // * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling
			cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'T'); // Bt*lambda
		}

		dual_solution_out.clear();
		dual_solution_out.resize(cluster.x_prim_cluster1.size());

		for (size_t d = 0; d < cluster.x_prim_cluster1.size(); d++) {
			dual_solution_out[d] = *cluster.x_prim_cluster1[d];
		}
		//dual_solution_out = cluster.x_prim_cluster1;

		for (size_t d = 0; d < cluster.domains.size(); d++) {

#ifdef BEM4I_TO_BE_REMOVED
			cluster.domains[d]->K.DenseMatVec(primal_solution_out[d], *cluster.x_prim_cluster2[d],'N');
#else
			cluster.domains[d]->K.MatVec(primal_solution_out[d], *cluster.x_prim_cluster2[d],'N');
#endif

			if (cluster.domains[d]->_RegMat.nnz > 0) {
				cluster.domains[d]->_RegMat.MatVecCOO(primal_solution_out[d], *cluster.x_prim_cluster2[d],'N', 1.0, -1.0); // K*u
			}

		}

		for (size_t d = 0; d < cluster.domains.size(); d++) {
			for (size_t pi = 0; pi < primal_solution_out[d].size(); pi++ ) {
				(*cluster.x_prim_cluster1[d])[pi] = in_right_hand_side_primal[d][pi]
										          - (*cluster.x_prim_cluster1[d])[pi];		// f - Bt*lambda


				KKT1_norm2[d] += (*cluster.x_prim_cluster1[d])[pi] * (*cluster.x_prim_cluster1[d])[pi]; // norm (f - Bt*lambda)

				(*cluster.x_prim_cluster1[d])[pi] = (*cluster.x_prim_cluster2[d])[pi] - (*cluster.x_prim_cluster1[d])[pi]; // K*u - (f - bt*lambda)
				KKT1_norm[d] += (*cluster.x_prim_cluster1[d])[pi] * (*cluster.x_prim_cluster1[d])[pi]; // norm (K*u - (f - bt*lambda))
			}
			KKT1_norm_cluster_local += KKT1_norm[d];
			KKT1_norm_cluster_local2 += KKT1_norm2[d];

		}

		MPI_Allreduce( &KKT1_norm_cluster_local, &KKT1_norm_cluster_global, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
		KKT1_norm_cluster_global = sqrt(KKT1_norm_cluster_global);

		MPI_Allreduce( &KKT1_norm_cluster_local2, &KKT1_norm_cluster_global2, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
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
			SEQ_VECTOR <double> tmp_dual(cluster.domains[d]->B1_comp_dom.rows, 0.0);
			cluster.domains[d]->B1_comp_dom.MatVec (primal_solution_out[d], tmp_dual, 'N', 0, 0, 0.0);
			for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
				cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += tmp_dual[i];
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

		MPI_Allreduce( &max_Bn_l, &max_Bn_l_g, 1, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);
		max_Bn_l_g = fabs(max_Bn_l_g);

		MPI_Allreduce( &lambda_n_max, &lambda_n_max_g, 1, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);
		if (fabs(lambda_n_max_g) < 10e-8)
			lambda_n_max_g += 1.0;


		MPI_Allreduce( &lambda_n_max_2, &lambda_n_max_2_g, 1, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);
		lambda_n_max_2_g = fabs(lambda_n_max_2_g);

		norm_ce  = parallel_norm_compressed(cluster, ce_l);
		if (fabs(norm_ce) < 1)
			norm_ce += 1.0;

		norm_cn  = parallel_norm_compressed(cluster, cn_l);
		if (fabs(norm_cn) < 1)
			norm_cn += 1.0;

		norm_Beu = parallel_norm_compressed(cluster, Be_l);
		parallel_norm_compressed(cluster, Bn_l);
		norm_Bn_lLambda = parallel_norm_compressed(cluster, Bn_lLambda);



		ESINFO(CONVERGENCE) << " **** Karush-Kuhn-Tucker-conditions ****     ";


		ESINFO(CONVERGENCE) << " Solution norm: norm(K*u - f + Bt*Lambda) =      				" << std::setw(6) << (KKT1_norm_cluster_global2 ? KKT1_norm_cluster_global / KKT1_norm_cluster_global2 :KKT1_norm_cluster_global);

		ESINFO(CONVERGENCE) << " Equality constraints: norm(Be*u - ce) =					" << std::setw(6) << norm_Beu / norm_ce;

		// KKT3
		ESINFO(CONVERGENCE) << " Inequality constraints: norm(max(Bn*u - cn,0)) =				" << std::setw(6) << max_Bn_l_g / norm_cn;
		ESINFO(CONVERGENCE) << " Check Multipliers positiveness: norm(max(-Lambda_N,0) =			" << std::setw(6) << lambda_n_max_2_g / lambda_n_max_g;
		ESINFO(CONVERGENCE) << " Check norm of Normal Multipliers: norm((Be*u - ce)*Lambda_N) =			" << std::setw(6) << norm_Bn_lLambda / ( norm_cn * lambda_n_max_g );

//		switch (configuration.solver) {
//		case ESPRESO_ITERATIVE_SOLVER::QPCE:
//			Solve_QPCE_singular_dom(cluster, in_right_hand_side_primal );
//			break;
//		default:
//			ESINFO(GLOBAL_ERROR) << "Unknown CG solver";
//		}


	}
}

void IterSolverBase::MakeSolution_Primal_singular_parallel (
		SuperCluster & cluster,
		SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal,
		SEQ_VECTOR < SEQ_VECTOR <double> > & primal_solution_out )  {

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_IterSolver_MakeSolution, "IterSolver_MakeSolution", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	primal_solution_out.clear();

	// R * mu
	SEQ_VECTOR<SEQ_VECTOR<double> > R_mu_prim_cluster;

	int amp_offset = 0;
	R_mu_prim_cluster.resize(cluster.domains.size());

	for (int c = 0; c < cluster.clusters.size(); c++) {
		for (size_t d = 0; d < cluster.clusters[c].domains.size(); d++) {
			SEQ_VECTOR <double > tmp (cluster.clusters[c].domains[d].domain_prim_size);

			if (USE_HFETI == 1) {
				if ( configuration.regularization == REGULARIZATION::FIX_POINTS ) {
					cluster.clusters[c].domains[d].Kplus_R .DenseMatVec(amplitudes, tmp, 'N', amp_offset, 0);
				} else {
					cluster.clusters[c].domains[d].Kplus_Rb.DenseMatVec(amplitudes, tmp, 'N', amp_offset, 0);
				}
			} else {
				cluster.clusters[c].domains[d].Kplus_R.DenseMatVec(     amplitudes, tmp, 'N', amp_offset, 0);
				amp_offset += cluster.clusters[c].domains[d].Kplus_R.cols;
			}

			R_mu_prim_cluster[ cluster.clusters[c].domains[d].domain_global_index ].swap(tmp);
		} // end domain loop
		if (USE_HFETI == 1) { amp_offset += cluster.clusters[c].G1.rows;}
	} // end cluster loop


	for (size_t d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d]->B1_comp_dom.rows );
		SEQ_VECTOR < double > tmp      ( cluster.domains[d]->domain_prim_size  );

		for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
			x_in_tmp[i] = dual_soultion_compressed_parallel[ cluster.domains[d]->lambda_map_sub_local[i]];
		}

		cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, tmp, 'T');

		for (size_t i = 0; i < tmp.size(); i++) {
			tmp[i] = in_right_hand_side_primal[d][i] - tmp[i];
		}
		primal_solution_out.push_back(tmp);
	}

	cluster.multKplus(primal_solution_out);

	#pragma omp parallel for
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		for (size_t i = 0; i < primal_solution_out[d].size()	; i++) {
			primal_solution_out[d][i] = primal_solution_out[d][i] + R_mu_prim_cluster[d][i];
		}
	}

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_IterSolver_MakeSolution);
#endif

}


// POWER Method
double IterSolverBase::Solve_power_method ( SuperCluster & cluster, double tol, eslocal maxit, eslocal method)
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
			Projector_Inv( timeEvalProj, cluster, Y, V , 0);
		} else {
			Projector    ( timeEvalProj, cluster, Y, V , 0);
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Y);


		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, Y, V, 0);
		} else {
			Projector    ( timeEvalProj, cluster, Y, V, 0);
		}

	}else{

					if (USE_GGtINV == 1) {
						Projector_Inv( timeEvalProj, cluster, Y, V , 0);
					} else {
						Projector    ( timeEvalProj, cluster, Y, V , 0);
					}

					apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Z);

					if (USE_GGtINV == 1) {
						Projector_Inv( timeEvalProj, cluster, Z, X, 0);
					} else {
						Projector    ( timeEvalProj, cluster, Z, X, 0);
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
					Projector_Inv( timeEvalProj, cluster, Y, V , 0);
				} else {
					Projector    ( timeEvalProj, cluster, Y, V , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Y);


				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, Y, V, 0);
				} else {
					Projector    ( timeEvalProj, cluster, Y, V, 0);
				}

    	}else{

				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, Y, V , 0);
				} else {
					Projector    ( timeEvalProj, cluster, Y, V , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, V, Z);

				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, Z, X, 0);
				} else {
					Projector    ( timeEvalProj, cluster, Z, X, 0);
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
//	MPI_Allreduce(&norm_x_l, &norm_x, 1, MPI_DOUBLE, MPI_MAX, environment->MPICommunicator);

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
void IterSolverBase::Solve_QPCE_singular_dom ( SuperCluster & cluster,
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
		Projector_Inv( timeEvalProj, cluster, b_l_, b_l , 0);
	} else {
		Projector    ( timeEvalProj, cluster, b_l_, b_l , 0);
	}
	// END*** projection of right hand side b



	// BEGIN*** Homogenization of the equality constraints and initialization
	if (USE_GGtINV == 1) {
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_im, 1 );
	} else {
		Projector	 ( timeEvalProj, cluster, cluster.vec_d, x_im, 1 );
	}

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_im, Ax_im);

	if (USE_GGtINV == 1) {
		Projector_Inv( timeEvalProj, cluster, Ax_im, tmp , 0);
	} else {
		Projector    ( timeEvalProj, cluster, Ax_im, tmp , 0);
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
			Projector_Inv( timeEvalProj, cluster, x_l, tmp , 0);
		} else {
			Projector    ( timeEvalProj, cluster, x_l, tmp , 0);
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
		} else {
			Projector    ( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
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
			case ESPRESO_PRECONDITIONER::LUMPED:
			case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
			case ESPRESO_PRECONDITIONER::DIRICHLET:
			case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
			case ESPRESO_PRECONDITIONER::MAGIC:

				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, p_l, w_l, 0 );
				} else {
					Projector( timeEvalProj, cluster, p_l, w_l, 0 );
				}

				Apply_Prec(timeEvalPrec, cluster, w_l, tmp_2);

				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, tmp_2, y_l, 0 );
				} else {
					Projector		  ( timeEvalProj, cluster, tmp_2, y_l, 0 );
				}

				for (size_t k = 0; k < p_l.size(); k++){
					p_l[k] = (y_l[k] * maxeig + 1.0 / rho * ( p_l[k]-w_l[k])) * _free[k];

				}

				break;
			case ESPRESO_PRECONDITIONER::NONE:


				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
			}








	cluster.G1_comp.MatVec(x_l, Cx_l, 'N');

	normx_l = parallel_norm_compressed(cluster, x_l);

	if (USE_GGtINV == 1) {
		Projector_Inv( timeEvalProj, cluster, Cx_l, tmp, 1 );
	} else {
		Projector	 ( timeEvalProj, cluster, Cx_l, tmp, 1 );
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
					Projector_Inv( timeEvalProj, cluster, p_l, tmp , 0);
				} else {
					Projector    ( timeEvalProj, cluster, p_l, tmp , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
				} else {
					Projector    ( timeEvalProj, cluster, Ax_l, PAPx_l, 0);
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

				MPI_Allreduce(&cnt_l, &cnt, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);

				if (cnt == 0) {
					alpha_f = INFINITY;
				} else {
					alpha_f_l = *std::min_element( tmp.begin(), tmp.begin() + cnt_l );
					MPI_Allreduce(&alpha_f_l, &alpha_f, 1, MPI_DOUBLE, MPI_MIN, environment->MPICommunicator);
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
							case ESPRESO_PRECONDITIONER::LUMPED:
							case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
							case ESPRESO_PRECONDITIONER::DIRICHLET:
							case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
							case ESPRESO_PRECONDITIONER::MAGIC:

								if (USE_GGtINV == 1) {
									Projector_Inv( timeEvalProj, cluster, tmp, w_l, 0 );
								} else {
									Projector( timeEvalProj, cluster, tmp, w_l, 0 );
								}

								Apply_Prec(timeEvalPrec, cluster, w_l, tmp_2);

								if (USE_GGtINV == 1) {
									Projector_Inv( timeEvalProj, cluster, tmp_2, y_l, 0 );
								} else {
									Projector		  ( timeEvalProj, cluster, tmp_2, y_l, 0 );
								}

								for (size_t k = 0; k < tmp.size(); k++){
									tmp[k] = (y_l[k] * maxeig + 1.0 / rho * ( tmp[k]-w_l[k])) * _free[k];

								}



								break;
							case ESPRESO_PRECONDITIONER::NONE:


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
							Projector_Inv( timeEvalProj, cluster, x_l, tmp , 0);
						} else {
							Projector    ( timeEvalProj, cluster, x_l, tmp , 0);
						}

						apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

						if (USE_GGtINV == 1) {
							Projector_Inv( timeEvalProj, cluster, Ax_l, g_l, 0);
						} else {
							Projector    ( timeEvalProj, cluster, Ax_l, g_l, 0);
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
							Projector_Inv( timeEvalProj, cluster, x_l, tmp , 0);
						} else {
							Projector    ( timeEvalProj, cluster, x_l, tmp , 0);
						}

						apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

						if (USE_GGtINV == 1) {
							Projector_Inv( timeEvalProj, cluster, Ax_l, g_l, 0);
						} else {
							Projector    ( timeEvalProj, cluster, Ax_l, g_l, 0);
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
					Projector_Inv( timeEvalProj, cluster, x_l, tmp , 0);
				} else {
					Projector    ( timeEvalProj, cluster, x_l, tmp , 0);
				}

				apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp, Ax_l);

				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, Ax_l, g_l, 0);
				} else {
					Projector    ( timeEvalProj, cluster, Ax_l, g_l, 0);
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
				Projector_Inv( timeEvalProj, cluster, Cx_l, tmp, 1 );
			} else {
				Projector	 ( timeEvalProj, cluster, Cx_l, tmp, 1 );
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
				Projector_Inv( timeEvalProj, cluster, mu, tmp, 1 );
			} else {
				Projector	 ( timeEvalProj, cluster, mu, tmp, 1 );
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
		Projector_Inv ( timeEvalProj, cluster, mu_tmp, amplitudes, 3 );
	} else {
		//TODO: Neni implementovan parametr 3
		Projector	  ( timeEvalProj, cluster, mu_tmp, amplitudes, 3 );
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
//	MPI_Allreduce(&cnt_l, &cnt, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);
//	MPI_Allreduce(&dnt_l, &dnt, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);


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



void IterSolverBase::Solve_RegCG ( SuperCluster & cluster,
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
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
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
		norm_prim_fl += cluster.domains[d]->norm_f;

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
	//MPI_Reduce   (&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, 0, environment->MPICommunicator);
	norm_prim_fg = sqrt(norm_prim_fg);



	// *** r = b - Ax *************************************************************
	#pragma omp parallel for
	for (size_t i = 0; i < r_l.size(); i++)
		r_l[i] = b_l[i] - Ax_l[i];

	if (USE_GGtINV == 1) {
		Projector_Inv( timeEvalProj, cluster, r_l, u_l , 0);
	} else {
		Projector    ( timeEvalProj, cluster, r_l, u_l , 0);
	}

	// *** Calculate the stop condition *******************************************
	//tol = epsilon * parallel_norm_compressed(cluster, u_l);

	double tol1 = epsilon * parallel_norm_compressed(cluster, u_l);
	double tol2 = epsilon * parallel_norm_compressed(cluster, b_l);

	if (tol1 < tol2 )
		tol = tol1;
	else
		tol = tol2;





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

//	double min_tol = 1e-12;
//	if (tol < min_tol) {
//		ESINFO(CONVERGENCE) << Info::TextColor::RED << "The NORM is fulfilled.";
//	} else {
		ESINFO(CONVERGENCE)
			<< spaces(indent.size() + iterationWidth - 4) << "iter"
			<< spaces(indent.size() + precision - 3) << "|r|" << spaces(2)
			<< spaces(indent.size() + 4) << "r" << spaces(4)
			<< spaces(indent.size() + (precision + 2) / 2 + (precision + 2) % 2 - 1) << "e" << spaces(precision / 2)
			<< spaces(indent.size()) << "time[s]";
//	}

	// *** Start the CG iteration loop ********************************************

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_RegCG_AllIterations, "RegCG-AllIterations", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	//for (int iter = 0; tol > min_tol && iter < CG_max_iter; iter++) {
	for (int iter = 0; iter < CG_max_iter; iter++) {

#ifdef READEX_LEVEL_1
		READEX_REGION_START(REG_RegCG_OneIteration, "RegCG-OneIteration", SCOREP_USER_REGION_TYPE_LOOP);
#endif

		timing.totalTime.start();

		#pragma omp parallel for
		for (size_t i = 0; i < r_l.size(); i++) {
			wp_l[i] = w_l[i];				//	wp = w;
			yp_l[i] = y_l[i];				//	yp = y
		}

		switch (USE_PREC) {
		case ESPRESO_PRECONDITIONER::LUMPED:
		case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
		case ESPRESO_PRECONDITIONER::DIRICHLET:
		case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
		case ESPRESO_PRECONDITIONER::MAGIC:
			proj1_time.start();
			if (USE_GGtINV == 1) {
				Projector_Inv( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector		  ( timeEvalProj, cluster, r_l, w_l, 0 );
			}
			proj1_time.end();

			// Scale
			prec_time.start();
			Apply_Prec(timeEvalPrec, cluster, w_l, z_l);
			prec_time.end();
			// Re-Scale

			proj2_time.start();
			if (USE_GGtINV == 1) {
				Projector_Inv( timeEvalProj, cluster, z_l, y_l, 0 );
			} else {
				Projector		  ( timeEvalProj, cluster, z_l, y_l, 0 );
			}
			proj2_time.end();
			break;
		case ESPRESO_PRECONDITIONER::NONE:
			proj_time.start();
			if (USE_GGtINV == 1) {
				Projector_Inv( timeEvalProj, cluster, r_l, w_l, 0 );
			} else {
				Projector		  ( timeEvalProj, cluster, r_l, w_l, 0 );
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

		//MPI_Allreduce(&norm_prim_l, &norm_prim_g, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
		////MPI_Reduce(&norm_prim_l, &norm_prim_g, 1, MPI_DOUBLE, MPI_SUM, 0, environment->MPICommunicator);
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
#ifdef READEX_LEVEL_1
		READEX_REGION_STOP(REG_RegCG_OneIteration);
#endif
		if (norm_l < tol)
			break;

	} // end of CG iterations

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_RegCG_AllIterations);
#endif

	// *** save solution - in dual and amplitudes *********************************************
	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;

	if (USE_GGtINV == 1) {
		Projector_Inv ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case ESPRESO_PRECONDITIONER::NONE:
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

void IterSolverBase::Solve_new_CG_singular_dom ( SuperCluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                            C G   -   N E W    I M P L E M E N T A T I O N                          # 
####################################################################################################*/
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
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d]->norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ax_l[i] - b_l[i];
  }

	if (USE_GGtINV == 1) {
		Projector_Inv( timeEvalProj, cluster, g_l, Pg_l , 0);
	} else {
		Projector    ( timeEvalProj, cluster, g_l, Pg_l , 0);
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
    case ESPRESO_PRECONDITIONER::LUMPED:
    case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
    case ESPRESO_PRECONDITIONER::DIRICHLET:
	  case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
    case ESPRESO_PRECONDITIONER::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, g_l, Pg_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, g_l, Pg_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      Apply_Prec(timeEvalPrec, cluster, Pg_l, MPg_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, MPg_l, z_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, MPg_l, z_l, 0 );
      }
      proj2_time.end();
      break;
    case ESPRESO_PRECONDITIONER::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, g_l, z_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, g_l, z_l, 0 );
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
		Projector_Inv ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case ESPRESO_PRECONDITIONER::NONE:
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

void IterSolverBase::Solve_full_ortho_CG_singular_dom ( SuperCluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                              C G      F U L L    O R T H O G O N A L                               # 
####################################################################################################*/
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
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d]->norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ax_l[i] - b_l[i];
  }

	if (USE_GGtINV == 1) {
		Projector_Inv( timeEvalProj, cluster, g_l, Pg_l , 0);
	} else {
		Projector    ( timeEvalProj, cluster, g_l, Pg_l , 0);
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
    case ESPRESO_PRECONDITIONER::LUMPED:
    case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
    case ESPRESO_PRECONDITIONER::DIRICHLET:
	  case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
    case ESPRESO_PRECONDITIONER::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, g_l, Pg_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, g_l, Pg_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      Apply_Prec(timeEvalPrec, cluster, Pg_l, MPg_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, MPg_l, z_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, MPg_l, z_l, 0 );
      }
      proj2_time.end();
      break;
    case ESPRESO_PRECONDITIONER::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, g_l, z_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, g_l, z_l, 0 );
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

	    MPI_Allreduce( &_Gamma_l[0], &Gamma_l[0], iter, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
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
  if (cnt_iter>0 && cnt_iter<1000 && cond_numb_FETI_operator && environment->MPIrank==0){
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
		Projector_Inv ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, g_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case ESPRESO_PRECONDITIONER::NONE:
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

void IterSolverBase::Solve_GMRES_singular_dom ( SuperCluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                                            G M R E S                                               # 
####################################################################################################*/
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
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d]->norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ax_l[i] - b_l[i];
  }

  switch (USE_PREC) {
  case ESPRESO_PRECONDITIONER::LUMPED:
  case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
  case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
  case ESPRESO_PRECONDITIONER::MAGIC:
    proj1_time.start();
    if (USE_GGtINV == 1) {
      Projector_Inv( timeEvalProj, cluster, g_l, Pg_l, 0 );
    } else {
      Projector		  ( timeEvalProj, cluster, g_l, Pg_l, 0 );
    }
    proj1_time.end();

    // Scale
    prec_time.start();
    Apply_Prec(timeEvalPrec, cluster, Pg_l, MPg_l);
    prec_time.end();
    // Re-Scale

    proj2_time.start();
    if (USE_GGtINV == 1) {
      Projector_Inv( timeEvalProj, cluster, MPg_l, z_l, 0 );
    } else {
      Projector		  ( timeEvalProj, cluster, MPg_l, z_l, 0 );
    }
    proj2_time.end();
    break;
  case ESPRESO_PRECONDITIONER::NONE:
    proj_time.start();
    if (USE_GGtINV == 1) {
      Projector_Inv( timeEvalProj, cluster, g_l, z_l, 0 );
    } else {
      Projector		  ( timeEvalProj, cluster, g_l, z_l, 0 );
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
    case ESPRESO_PRECONDITIONER::LUMPED:
    case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
    case ESPRESO_PRECONDITIONER::DIRICHLET:
	  case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
    case ESPRESO_PRECONDITIONER::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, w_l, Pw_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, w_l, Pw_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      Apply_Prec(timeEvalPrec, cluster, Pw_l, MPw_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, MPw_l, z_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, MPw_l, z_l, 0 );
      }
      proj2_time.end();
      break;
    case ESPRESO_PRECONDITIONER::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, w_l, z_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, w_l, z_l, 0 );
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
		Projector_Inv ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, w_l, amplitudes, 2 );
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
		Projector_Inv ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case ESPRESO_PRECONDITIONER::NONE:
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
void IterSolverBase::Solve_BICGSTAB_singular_dom ( SuperCluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                                         B I C G S T A B                                            # 
####################################################################################################*/
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
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	} else {
		Projector	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	}

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ay_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	double norm_prim_fl = 0.0;
	double norm_prim_fg = 0.0;
	for (size_t d = 0; d < cluster.domains.size(); d++){
		norm_prim_fl += cluster.domains[d]->norm_f;
  }

	MPI_Allreduce(&norm_prim_fl, &norm_prim_fg, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
	norm_prim_fg = sqrt(norm_prim_fg);

	// *** g = Ax - b *************************************************************
	#pragma omp parallel for
for (size_t i = 0; i < g_l.size(); i++){
		g_l[i] = Ay_l[i] - b_l[i];
  }

  switch (USE_PREC) {
  case ESPRESO_PRECONDITIONER::LUMPED:
  case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
  case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
  case ESPRESO_PRECONDITIONER::MAGIC:
    proj1_time.start();
    if (USE_GGtINV == 1) {
      Projector_Inv( timeEvalProj, cluster, g_l, tmp1_l, 0 );
    } else {
      Projector		  ( timeEvalProj, cluster, g_l, tmp1_l, 0 );
    }
    proj1_time.end();

    // Scale
    prec_time.start();
    Apply_Prec(timeEvalPrec, cluster, tmp1_l, g_l);
    prec_time.end();
    // Re-Scale

    proj2_time.start();
    if (USE_GGtINV == 1) {
      Projector_Inv( timeEvalProj, cluster, g_l, z_l, 0 );
    } else {
      Projector		  ( timeEvalProj, cluster, g_l, z_l, 0 );
    }
    proj2_time.end();
    break;
  case ESPRESO_PRECONDITIONER::NONE:
    proj_time.start();
    if (USE_GGtINV == 1) {
      Projector_Inv( timeEvalProj, cluster, g_l, z_l, 0 );
    } else {
      Projector		  ( timeEvalProj, cluster, g_l, z_l, 0 );
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
    case ESPRESO_PRECONDITIONER::LUMPED:
    case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
    case ESPRESO_PRECONDITIONER::DIRICHLET:
	  case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
    case ESPRESO_PRECONDITIONER::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, Ay_l, v_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, Ay_l, v_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      Apply_Prec(timeEvalPrec, cluster, v_l, tmp1_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, tmp1_l, v_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, tmp1_l, v_l, 0 );
      }
      proj2_time.end();
      break;
    case ESPRESO_PRECONDITIONER::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, Ay_l, v_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, Ay_l, v_l, 0 );
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
    case ESPRESO_PRECONDITIONER::LUMPED:
    case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
    case ESPRESO_PRECONDITIONER::DIRICHLET:
	  case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
    case ESPRESO_PRECONDITIONER::MAGIC:
      proj1_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, Ay_l, t_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, Ay_l, t_l, 0 );
      }
      proj1_time.end();

      // Scale
      prec_time.start();
      Apply_Prec(timeEvalPrec, cluster, t_l, tmp1_l);
      prec_time.end();
      // Re-Scale

      proj2_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, tmp1_l, t_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, tmp1_l, t_l, 0 );
      }
      proj2_time.end();
      break;
    case ESPRESO_PRECONDITIONER::NONE:
      proj_time.start();
      if (USE_GGtINV == 1) {
        Projector_Inv( timeEvalProj, cluster, Ay_l, t_l, 0 );
      } else {
        Projector		  ( timeEvalProj, cluster, Ay_l, t_l, 0 );
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
		Projector_Inv ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj , cluster, w_l, amplitudes, 2 );
	}


  dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = w_l;


	if (USE_GGtINV == 1) {
		Projector_Inv ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, w_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case ESPRESO_PRECONDITIONER::NONE:
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



void IterSolverBase::Solve_PipeCG_singular_dom ( SuperCluster & cluster,
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
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );
	else
		Projector	 ( timeEvalProj, cluster, cluster.vec_d, x_l, 1 );


	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);


	// *** Ax = apply_A(CLUSTER,Bt,x); ********************************************
	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l); //apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);


	// *** r = b - Ax; ************************************************************
	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:

		#pragma omp parallel for
for (size_t i = 0; i < r_l.size(); i++) {
			tmp_l[i] = b_l[i] - Ax_l[i];
		}

		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, tmp_l, r_l, 0 );
		} else {
			Projector    ( timeEvalProj, cluster, tmp_l, r_l, 0 );
		}

		tol = epsilon * parallel_norm_compressed(cluster, r_l);

		Apply_Prec(timeEvalPrec, cluster, r_l, tmp_l);
		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, tmp_l, u_l, 0 );
		} else {
			Projector    ( timeEvalProj, cluster, tmp_l, u_l, 0 );
		}

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, u_l, tmp_l); //apply_A_l_compB(timeEvalAppa, cluster, u_l, tmp_l);
		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, tmp_l, w_l, 0 );
		} else {
			Projector    ( timeEvalProj, cluster, tmp_l, w_l, 0 );
		}

		break;
	case ESPRESO_PRECONDITIONER::NONE:
		#pragma omp parallel for
		for (size_t i = 0; i < r_l.size(); i++) {
			r_l[i] = b_l[i] - Ax_l[i];
		}

		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, r_l, u_l, 0 );
		} else {
			Projector    ( timeEvalProj, cluster, r_l, u_l, 0 );
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
		case ESPRESO_PRECONDITIONER::LUMPED:
		case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
		case ESPRESO_PRECONDITIONER::DIRICHLET:
	  case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
		case ESPRESO_PRECONDITIONER::MAGIC:
			parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, r_l, &mpi_req, reduction_tmp, send_buf); // norm_l = parallel_norm_compressed(cluster, r_l);
			ddot_time.end();

			prec_time.start();
			Apply_Prec(timeEvalPrec, cluster, w_l, tmp_l);
			prec_time.end();

			proj_time.start();
			if (USE_GGtINV == 1) {
				Projector_Inv( timeEvalProj, cluster, tmp_l, m_l, 0 );
			} else {
				Projector    ( timeEvalProj, cluster, tmp_l, m_l, 0 );
			}
			proj_time.end();

			appA_time.start();
			apply_A_l_comp_dom_B(timeEvalAppa, cluster, m_l, tmp_l); //apply_A_l_compB(timeEvalAppa, cluster, m_l, tmp_l);
			appA_time.end();

			proj_time.start();
			if (USE_GGtINV == 1) {
				Projector_Inv( timeEvalProj, cluster, tmp_l, n_l, 0 );
			} else {
				Projector    ( timeEvalProj, cluster, tmp_l, n_l, 0 );
			}
			proj_time.end();

			break;
		case ESPRESO_PRECONDITIONER::NONE:
			parallel_ddot_compressed_non_blocking(cluster, r_l, u_l, w_l, u_l, u_l, &mpi_req, reduction_tmp, send_buf); // norm_l = parallel_norm_compressed(cluster, u_l);
			ddot_time.end();

			proj_time.start();

			if (USE_GGtINV == 1) {
				Projector_Inv( timeEvalProj, cluster, w_l, m_l, 0 );
			} else {
				Projector    ( timeEvalProj, cluster, w_l, m_l, 0 );
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
		Projector_Inv ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************
	timing.addEvent(ddot_time);

	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:
		//timing.addEvent(proj1_time);
		timing.addEvent(proj_time);
		timing.addEvent(prec_time );
		//timing.addEvent(proj2_time);
		break;
	case ESPRESO_PRECONDITIONER::NONE:
		timing.addEvent(proj_time);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
	}

	timing.addEvent(appA_time);
	timing.addEvent(vec_time );

}

void IterSolverBase::CreateConjProjector(Cluster & cluster) {
/*

	//int d = 0;


	cluster.SetupPreconditioner();

	SEQ_VECTOR<SparseMatrix> Ml_per_dom (cluster.domains.size());
	SEQ_VECTOR<SparseMatrix> Ml_fin_pd  (cluster.domains.size());

	SEQ_VECTOR<SparseMatrix> Bdec   (cluster.domains.size());
	SEQ_VECTOR<SparseMatrix> Btdec  (cluster.domains.size());



	for (size_t d = 0; d < cluster.domains.size(); d++) {


		cluster.domains[d].B1 = cluster.instance->B1[d];

		SparseMatrix B1full = cluster.instance->B1[d];
		B1full.ConvertToCSR(0);

		SparseMatrix B1tfull = cluster.instance->B1[d];
		B1tfull.ConvertToCSR(0);
		B1tfull.MatTranspose();

		SparseMatrix PrecFull = cluster.domains[d].Prec;

		PrecFull.ConvertDenseToCSR(1);


		SparseMatrix Bt = cluster.domains[d].B1t_DirPr;
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
		prec_map_vec = cluster.domains[d].B1t_Dir_perm_vec;

		for (size_t i = 0; i < PrecFull.I_row_indices.size(); i++) {
			PrecFull.I_row_indices[i] = prec_map_vec[PrecFull.I_row_indices[i]-1];
			PrecFull.J_col_indices[i] = prec_map_vec[PrecFull.J_col_indices[i]-1];
		}

		PrecFull.rows = cluster.domains[d].B1.cols;
		PrecFull.cols = cluster.domains[d].B1.cols;
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

		std::cout << ".";

	}

	std::cout << std::endl;

	//TODO: can be done as log N
	SparseMatrix Ml_per_cluster;
	Ml_per_cluster = Ml_per_dom[0];
	for (size_t d = 1; d < cluster.domains.size(); d++) {
		Ml_per_cluster.MatAddInPlace(Ml_per_dom[d],'N',1.0);
	}

	SparseMatrix C_per_cluster;

	//std::cout << Ml_per_cluster.SpyText();

	for (size_t d = 0; d < cluster.domains.size(); d++) {

		SparseMatrix B1tfull = cluster.instance->B1[d];
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
		//std::cout << Ml_fin_pd[d].SpyText();
		Ml_fin_pd[d].ConvertToCOO(0);

		SparseMatrix Prec = cluster.domains[d].Prec;
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



		int il = cluster.domains[d].Kplus_R.cols; 								// first eigen value to be calculated
		int iv = cluster.domains[d].Kplus_R.cols + configuration.GENEO_SIZE;	// last  eigen value to be calculated

		std::cout << "FETI Geneo eigen vectors from " << il << " to " << iv << std::endl;


		int eig_n = (iv-il+1);
		eslocal n = A.rows;
		SEQ_VECTOR <double> w (eig_n);				// eigen values storage
		SEQ_VECTOR <double> eig_vectors (eig_n*n);	// eigen vectors storage
		SEQ_VECTOR <eslocal> ifail (n);				// dummy
		SEQ_VECTOR <eslocal> m (n);					// dummy
		double tmpd;								// dummy
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


			std::cout << ".";

	}

	std::cout << std::endl;

	std::cout << "1" << std::endl;

	C_per_cluster.MatTranspose();
	C_per_cluster.ConvertCSRToDense(0);

	SparseMatrix CP_per_cluster;
	SparseMatrix CPt_per_cluster;
	SparseMatrix ACP_per_cluster;

	std::cout << "2" << std::endl;

	Projector_l_inv_compG(timeEvalProj, cluster, C_per_cluster, CP_per_cluster );

	std::cout << "3" << std::endl;

	apply_A_l_Mat (timeEvalAppa, cluster, CP_per_cluster, ACP_per_cluster);

	SparseMatrix Proj_per_cluster;
	SparseMatrix Proj_per_cluster2;

	SparseMatrix Proj_per_cluster_sp;

	CPt_per_cluster = CP_per_cluster;
	CP_per_cluster.ConvertDenseToCSR(0);

	CPt_per_cluster.ConvertDenseToCSR(1);
	CPt_per_cluster.MatTranspose();
	CPt_per_cluster.ConvertCSRToDense(1);

	std::cout << "4" << std::endl;

	Proj_per_cluster.DenseMatMat(CPt_per_cluster,'N',ACP_per_cluster,'N');

	Proj_per_cluster2.DenseMatMat(CP_per_cluster,'T',ACP_per_cluster,'N'); // TODO : DGEMM nefunguje s transpose turned on - needs to be fixed

	std::cout << "5" << std::endl;

	ACP_per_cluster.ConvertDenseToCSR(0);
	CPt_per_cluster.ConvertDenseToCSR(0);

	Proj_per_cluster_sp.MatMat(CPt_per_cluster,'N',ACP_per_cluster);
	Proj_per_cluster_sp.ConvertCSRToDense(0);


	Proj_per_cluster.RemoveLowerDense();
	Proj_per_cluster.mtype = espreso::MatrixType::REAL_SYMMETRIC_POSITIVE_DEFINITE;

	std::cout << "6" << std::endl;

	DenseSolverMKL ProjSolver;

	ProjSolver.ImportMatrix(Proj_per_cluster);

	ProjSolver.Factorization("Geneo projector factorization");

	std::cout << "7" << std::endl;

	// *** Pass objects to cluster
	cluster.CFCt.ImportMatrix(Proj_per_cluster);
	cluster.CFCt.Factorization("Geneo projector factorization");
	cluster.C   = CP_per_cluster;
	cluster.Ct  = CPt_per_cluster;


*/
}


void IterSolverBase::ConjProj(  Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {

	/*
	SEQ_VECTOR <double> tmp_x_in, tmp1, tmp2, tmp3, tmp4;

	tmp_x_in = x_in;

	tmp1.resize(x_in.size(), 0.0);
	tmp2.resize(x_in.size(), 0.0);
	tmp3.resize(x_in.size(), 0.0);
	tmp4.resize(x_in.size(), 0.0);

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp_x_in, tmp1);

	cluster.Ct.DenseMatVec(tmp1,tmp2);

	cluster.CFCt.Solve(tmp2,tmp3,1);

	cluster.C.DenseMatVec(tmp3, tmp4);

	for (size_t i = 0; i < x_in.size(); i++)
		y_out[i] = x_in[i] - tmp4[i];
*/
}


void IterSolverBase::ConjProj_t(Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {

/*
	SEQ_VECTOR <double> tmp_x_in, tmp1, tmp2, tmp3, tmp4;

	tmp_x_in = x_in;

	tmp1.resize(x_in.size(), 0.0);
	tmp2.resize(x_in.size(), 0.0);
	tmp3.resize(x_in.size(), 0.0);
	tmp4.resize(x_in.size(), 0.0);

	cluster.Ct.DenseMatVec(tmp_x_in,tmp1);

	cluster.CFCt.Solve(tmp1,tmp2,1);

	cluster.C.DenseMatVec(tmp2, tmp3);

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp3, tmp4);

	for (size_t i = 0; i < x_in.size(); i++)
		y_out[i] = x_in[i] - tmp4[i];
*/
}



void IterSolverBase::ConjProj_lambda0(  Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {
/*
	SEQ_VECTOR <double> tmp_x_in, tmp1, tmp2;

	tmp1.resize(x_in.size(), 0.0);
	tmp2.resize(x_in.size(), 0.0);

	cluster.Ct.  DenseMatVec(x_in, tmp1);
	cluster.CFCt.Solve		(tmp1, tmp2, 	1);
	cluster.C.   DenseMatVec(tmp2, y_out);
*/
}



void IterSolverBase::Solve_full_ortho_CG_singular_dom_geneo ( SuperCluster & cluster,
	    SEQ_VECTOR < SEQ_VECTOR <double> > & in_right_hand_side_primal)
{
/*####################################################################################################
#                              C G      F U L L    O R T H O G O N A L                               #
####################################################################################################*/
//


	std::cout << "\n\n Full orthogonal with restart and conjugate projector \n";

	if (configuration.conj_projector == CONJ_PROJECTOR::GENEO) {
		std::cout << " - using GENEO \n";
	}

	size_t dl_size = cluster.my_lamdas_indices.size();

	SEQ_VECTOR <double> x_l (dl_size, 0);
	SEQ_VECTOR <double> x_0 (dl_size, 0);
	SEQ_VECTOR <double> x_1 (dl_size, 0);

	SEQ_VECTOR <double> Ax_l(dl_size, 0);
	SEQ_VECTOR <double> Ax_0(dl_size, 0);


	SEQ_VECTOR <double> r_l(dl_size, 0);
	SEQ_VECTOR <double> r_0(dl_size, 0);
	SEQ_VECTOR <double> r_tmp(dl_size, 0);


	SEQ_VECTOR <double> Pg_l(dl_size, 0);
	SEQ_VECTOR <double> MPg_l(dl_size, 0);
	SEQ_VECTOR <double> z_l(dl_size, 0);
	SEQ_VECTOR <double> _z_l(dl_size, 0);
	SEQ_VECTOR <double> w_l(dl_size, 0);
	SEQ_VECTOR <double> _w_l(dl_size, 0);
	SEQ_VECTOR <double> Aw_l(dl_size, 0);
	SEQ_VECTOR <double> PAw_l(dl_size, 0);
	SEQ_VECTOR <double> b_l(dl_size, 0);
	SEQ_VECTOR <double> v_tmp_l(dl_size, 0);

//	SEQ_VECTOR <double> d_H(CG_max_iter, 0);
//	SEQ_VECTOR <double> e_H(CG_max_iter, 0);

	double rho_l;
	double rho_l_prew = 1;
	double norm_l = 1e33;
	double norm_l_new;
	bool norm_decrease = true;
	double tol;
	double ztg;
	double wtAw;
	int N_ITER=0 ;


	double delta, gamma, fi;
	SEQ_VECTOR <double> delta_l	 	(CG_max_iter, 0);
//	SEQ_VECTOR <double> gamma_l	 	(CG_max_iter, 0);
	SEQ_VECTOR <double> fi_l	 	(CG_max_iter, 0);
	SEQ_VECTOR <double> fi_g	 	(CG_max_iter, 0);

	bool update = true;

	int cnt_iter=0;

	int restart_num = 0;

	eslocal restart_iter = 0;


	cluster.CreateVec_b_perCluster ( in_right_hand_side_primal );	// prava strana dualu
	cluster.CreateVec_d_perCluster ( in_right_hand_side_primal );	// e - vektor tuhych pohybu

	// *** Combine vectors b from all clusters ************************************
	All_Reduce_lambdas_compB(cluster, cluster.vec_b_compressed, b_l);

	if (USE_GGtINV == 1) {
		Projector_Inv( timeEvalProj, cluster, cluster.vec_d, x_0, 1 );
	} else {
		Projector( timeEvalProj, cluster, cluster.vec_d, x_0, 1 );
	}

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_0, Ax_0);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	#pragma omp parallel for
	for (size_t i = 0; i < r_0.size(); i++){
		r_0[i] =  b_l[i] - Ax_0[i];
	}

//	if (configuration.conj_projector == CONJ_PROJECTOR::GENEO) {
//		ConjProj_lambda0(cluster, r_0, x_1);
//		ConjProj_t(cluster, r_0, r_0);
//	}

	SparseMatrix W_l;
	W_l.type = 'G';
	W_l.rows = dl_size;
	W_l.cols = 0;


	SparseMatrix AW_l;
	AW_l.type = 'G';
	AW_l.rows = dl_size;
	AW_l.cols = 0;

	SEQ_VECTOR <double> Gamma_l  	(CG_max_iter, 0);
	SEQ_VECTOR <double> _Gamma_l 	(CG_max_iter, 0);
	SEQ_VECTOR <double> WtAW_l	 	(CG_max_iter, 0);

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::DIRICHLET:

		proj1_time.start();
		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, r_0, r_l, 0 );
		} else {
			Projector		  ( timeEvalProj, cluster, r_0, r_l, 0 );
		}
		proj1_time.end();

		// Scale

		prec_time.start();
		Apply_Prec(timeEvalPrec, cluster, r_l, z_l);
		prec_time.end();

		// Re-Scale

		proj2_time.start();
		if (USE_GGtINV == 1) {
			Projector_Inv( timeEvalProj, cluster, z_l, w_l, 0 );
		} else {
			Projector		  ( timeEvalProj, cluster, z_l, w_l, 0 );
		}
		proj2_time.end();

		break;

	default:
		ESINFO(GLOBAL_ERROR) << "FETI Geneo requires dirichlet preconditioner.";
	}



	// *** Calculate the stop condition *******************************************

	double tol1;
//	if (configuration.conj_projector == CONJ_PROJECTOR::GENEO) {
//		ConjProj(cluster, r_l, r_tmp );
//		tol1 = epsilon * parallel_norm_compressed(cluster, r_tmp);
//	}
//	else{
		tol1 = epsilon * parallel_norm_compressed(cluster, r_l);
//	}

	double tol2 = epsilon * parallel_norm_compressed(cluster, b_l);

	if (tol1 < tol2 )
		tol = tol1;
	else
		tol = tol2;

	int precision      = ceil(log(1 / epsilon) / log(10)) + 1;
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

	// *** END - Calculate the stop condition *******************************************



	W_l.dense_values.insert(W_l.dense_values.end(), w_l.begin(), w_l.end());
	W_l.nnz += w_l.size();
	W_l.cols++;

	eslocal iter = 0;

	// *** Start the CG iteration loop ********************************************
	for (int t_iter = 0; t_iter < CG_max_iter; t_iter++) {

		timing.totalTime.start();

		N_ITER +=1;

		appA_time.start();
		apply_A_l_comp_dom_B(timeEvalAppa, cluster, w_l, Aw_l);
		appA_time.end();

//		if (configuration.conj_projector == CONJ_PROJECTOR::GENEO) {
//			ConjProj_t(cluster, Aw_l, Aw_l);
//		}

		AW_l.dense_values.insert(AW_l.dense_values.end(), Aw_l.begin(), Aw_l.end());
		AW_l.nnz += Aw_l.size();
		AW_l.cols++;

		delta = parallel_ddot_compressed(cluster, w_l, Aw_l);
		delta_l[iter] = delta;
		gamma = parallel_ddot_compressed(cluster, z_l, r_l);


		Projector_Inv(timeEvalProj, cluster, Aw_l, PAw_l, 0);

		#pragma omp parallel for
		for (size_t i = 0; i < x_l.size(); i++) {
			x_l[i] = x_l[i] + gamma/delta_l[iter] *   w_l[i];

            r_l[i] = r_l[i] - gamma/delta_l[iter] * PAw_l[i];

		}


		switch (USE_PREC) {
		case ESPRESO_PRECONDITIONER::DIRICHLET:

			// Scale

			prec_time.start();
			Apply_Prec(timeEvalPrec, cluster, r_l, z_l);
			prec_time.end();

			// Re-Scale

			proj2_time.start();
			if (USE_GGtINV == 1) {
				Projector_Inv( timeEvalProj, cluster, z_l, w_l, 0 );
			} else {
				Projector		  ( timeEvalProj, cluster, z_l, w_l, 0 );
			}
			proj2_time.end();

			break;

		default:
			ESINFO(GLOBAL_ERROR) << "FETI Geneo requires dirichlet preconditioner.";
		}


		// filtering duplicit Lambda entries
		#pragma omp parallel for
		for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++) {
			_w_l[i] = w_l[i] * cluster.my_lamdas_ddot_filter[i];
		}
		AW_l.DenseMatVec(_w_l,fi_l ,'T');

//		std::cout << std::endl;
//		for (int i = 0; i < AW_l.cols; i++) {
//			std::cout << fabs (fi_l[i]) << ", ";
//		}
//		std::cout << std::endl;

		update = true;
		if ( iter > configuration.RESTART_ITER) {

			double fi_last = fabs( fi_l[iter - 1] );
			if (restart_num < configuration.NUM_RESTARTS){
				for (eslocal i = 0; i < iter - 1; i++) {

					if (!norm_decrease){
						if ( fabs (fi_l[i]) < fi_last ) {
							update = false;
						}
					}
				}
		    }
		}


		if (update) {

			#pragma omp parallel for
			for (eslocal i = 0; i <= iter; i++) {
				fi_g[i] = fi_l[i] / delta_l[i];
			}

			MPI_Allreduce( &fi_g[0], &fi_l[0], iter + 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
			W_l.DenseMatVec(fi_l,v_tmp_l);

			#pragma omp parallel for
			for (size_t i = 0; i < x_l.size(); i++) {
				w_l[i] = w_l[i] -  v_tmp_l[i];
			}

			W_l.dense_values.insert(W_l.dense_values.end(), w_l.begin(), w_l.end());
			W_l.nnz += w_l.size();
			W_l.cols++;


			iter++;

		} else {

			iter = 0;

			restart_num += 1;


			ESINFO(CONVERGENCE) << "Restarting orthogonalization ... ";

			// Restart W_l
			W_l.dense_values.clear();
			W_l.nnz  = 0;
			W_l.cols = 0;

			// Restart AW_l
			AW_l.dense_values.clear();
			AW_l.nnz  = 0;
			AW_l.cols = 0;

			#pragma omp parallel for
			for (size_t i = 0; i < r_0.size(); i++){
				x_0[i] =  x_0[i] + x_l[i];
				x_l[i] = 0.0;
			}

			apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_0, Ax_0);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

			#pragma omp parallel for
			for (size_t i = 0; i < r_0.size(); i++){
				r_0[i] =  b_l[i] - Ax_0[i];
			}

//			if (configuration.conj_projector == CONJ_PROJECTOR::GENEO) {
//				ConjProj_lambda0(cluster, r_0, x_l);
//				ConjProj_t(cluster, r_0, r_0);
//			}

			switch (USE_PREC) {
			case ESPRESO_PRECONDITIONER::DIRICHLET:

				proj1_time.start();
				if (USE_GGtINV == 1) {
					Projector_Inv(timeEvalProj, cluster, r_0, r_l, 0 );
				} else {
					Projector		  ( timeEvalProj, cluster, r_0, r_l, 0 );
				}
				proj1_time.end();

				// Scale

				prec_time.start();
				Apply_Prec(timeEvalPrec, cluster, r_l, z_l);
				prec_time.end();

				// Re-Scale

				proj2_time.start();
				if (USE_GGtINV == 1) {
					Projector_Inv( timeEvalProj, cluster, z_l, w_l, 0 );
				} else {
					Projector		  ( timeEvalProj, cluster, z_l, w_l, 0 );
				}
				proj2_time.end();

				break;

			default:
				ESINFO(GLOBAL_ERROR) << "FETI Geneo requires dirichlet preconditioner.";
			}

			W_l.dense_values.insert(W_l.dense_values.end(), w_l.begin(), w_l.end());
			W_l.nnz += w_l.size();
			W_l.cols++;

		}

		norm_time.start();
		norm_l_new = sqrt (parallel_ddot_compressed(cluster, r_l, z_l) );

		if (iter > 2) {
			if (norm_l > norm_l_new){
				norm_decrease = true;
			}else{
				norm_decrease = false;
			}
		}

		norm_l = norm_l_new;
		norm_time.end();

		timing.totalTime.end();

		ESINFO(CONVERGENCE)
		<< indent << std::setw(iterationWidth) << iter
		<< indent << std::fixed << std::setprecision(precision) 	<< norm_l / tol * epsilon
		<< indent << std::scientific << std::setprecision(3) 		<< norm_l
		<< indent << std::fixed << std::setprecision(precision - 1) << epsilon
		<< indent << std::fixed << std::setprecision(5) 			<< timing.totalTime.getLastStat();

		// *** Stop condition ******************************************************************
		if (norm_l < tol)
			break;

	} // end of CG iterations

	ESINFO(CONVERGENCE) << "FULL CG Stop after " << N_ITER << " Ïterations";
	// *** save solution - in dual and amplitudes *********************************************

//	if (configuration.conj_projector == CONJ_PROJECTOR::GENEO) {
//		ConjProj(cluster, x_l, x_l);
//	}

	#pragma omp parallel for
	for (size_t i = 0; i < x_l.size(); i++) {
		x_l[i] = x_l[i] + x_0[i] + x_1[i];
	}

	apply_A_l_comp_dom_B(timeEvalAppa, cluster, x_l, Ax_l);// apply_A_l_compB(timeEvalAppa, cluster, x_l, Ax_l);

	#pragma omp parallel for
	for (size_t i = 0; i < x_l.size(); i++){
		r_l[i] =  b_l[i] - Ax_l[i];
	}


	dual_soultion_compressed_parallel   = x_l;
	dual_residuum_compressed_parallel   = r_l;


	if (USE_GGtINV == 1) {
		Projector_Inv ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	} else {
		Projector	  ( timeEvalProj, cluster, r_l, amplitudes, 2 );
	}
	// *** end - save solution - in dual and amplitudes ***************************************


	// *** Preslocal out the timing for the iteration loop ***************************************

	switch (USE_PREC) {
	case ESPRESO_PRECONDITIONER::LUMPED:
	case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
	case ESPRESO_PRECONDITIONER::DIRICHLET:
	case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
	case ESPRESO_PRECONDITIONER::MAGIC:
		timing.addEvent(proj1_time);
		timing.addEvent(prec_time );
		timing.addEvent(proj2_time);
		break;
	case ESPRESO_PRECONDITIONER::NONE:
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



















// *** Coarse problem routines *******************************************
void IterSolverBase::CreateGGt( SuperCluster & cluster )

{

	ESINFO(GLOBAL_ERROR) << "Projector with factorized GGt matrix is not supported in current version";

	// TODO: Obsolete code - must be updated before used with current version
    // Code is commented

//	SparseMatrix G;
//
//
//	//if (mpi_rank == mpi_root)
//	//	G.MatAppend(cluster.G1);
//
//	//for (eslocal mr = 1; mr < mpi_size; mr++) {
//	//	SparseMatrix Gtmp;
//	//	SendMatrix(mpi_rank, mr, cluster.G1, mpi_root, Gtmp);
//
//	//	if (mpi_rank == mpi_root) {
//	//		G.MatAppend(Gtmp);
//	//		Gtmp.Clear();
//	//	}
//	//}
//
//	//// **** Log N MPI reduce
//	eslocal count_cv = 0;
//	for (eslocal li = 2; li <= 2*mpi_size; li = li * 2 ) {
//
//		SparseMatrix recv_m;
//
//		if (mpi_rank % li == 0) {
//			if (li == 2)
//				G.MatAppend(cluster.G1);
//			if ((mpi_rank + li/2) < mpi_size) {
//				SendMatrix(mpi_rank, mpi_rank + li/2, cluster.G1, mpi_rank,        recv_m);
//				G.MatAppend(recv_m);
//			} else {
//				SendMatrix(mpi_rank, mpi_size + 1, cluster.G1, mpi_size + 1,        recv_m);
//			}
//		} else {
//
//			if ((mpi_rank + li/2) % li == 0)
//			{
//				if (li == 2)
//					SendMatrix(mpi_rank, mpi_rank       , cluster.G1, mpi_rank - li/2, recv_m);
//				else
//					SendMatrix(mpi_rank, mpi_rank       , G         , mpi_rank - li/2, recv_m);
//			} else {
//				SendMatrix(mpi_rank, mpi_rank+1, cluster.G1, mpi_rank+1,recv_m);
//			}
//		}
//
//		MPI_Barrier(environment->MPICommunicator);
//
//		count_cv += mpi_size/li;
//
//		ESINFO(PROGRESS3) << " Collecting matrices G : " << count_cv <<" of " << mpi_size;
//	}
//
//	//SparseMatrix Gtt;
//	if (mpi_rank != mpi_root)
//		G.Clear();
//	//else {
//	//Gtt = G;
//	//G.Clear();
//	//}
//	// ****
//
//	if (mpi_rank == mpi_root) {
//
//		MKL_Set_Num_Threads(PAR_NUM_THREADS);
//		// Create Gt and later GGt matrices and remove all elements under main diagonal of the GGt
//		SparseMatrix Gt;
//
//		double t1 = omp_get_wtime();
//		G.MatTranspose(Gt);
//		ESINFO(PROGRESS3) << "Gtranspose = " << omp_get_wtime() - t1;
//
//		t1 = omp_get_wtime();
//		SparseMatrix GGt_Mat;
//		GGt_Mat.MatMat(G, 'N', Gt);
//		ESINFO(PROGRESS3) << "G x Gt = " << omp_get_wtime() - t1;
//
//		t1 = omp_get_wtime();
//		Gt.Clear();
//		G.Clear();
//		ESINFO(PROGRESS3) << "G and Gt clear = " << omp_get_wtime() - t1;
//
//		ESINFO(EXHAUSTIVE) << GGt_Mat.SpyText();
//
//		t1 = omp_get_wtime();
//		GGt_Mat.RemoveLower();
//		ESINFO(PROGRESS3) << "GGt remove lower = " << omp_get_wtime() - t1;
//
//		t1 = omp_get_wtime();
//		// Create Sparse Direct solver for GGt
//		GGt.msglvl = Info::report(LIBRARIES) ? 1 : 0;
//
//		t1 = omp_get_wtime();
//		GGt.ImportMatrix(GGt_Mat);
//		ESINFO(PROGRESS3) << "ImportMatrix = " << omp_get_wtime() - t1;
//
//
//		t1 = omp_get_wtime();
//		GGt_Mat.Clear();
//
//
//		t1 = omp_get_wtime();
//		std::stringstream ss;
//		ss << "Create GGt -> rank: " << environment->MPIrank;
//		GGt.Factorization(ss.str());
//		ESINFO(PROGRESS3) << "Factorization = " << omp_get_wtime() - t1;
//
//
//		t1 = omp_get_wtime();
//		GGt.msglvl = 0;
//		//TODO:
//		MKL_Set_Num_Threads(1);
//	}
//
//
//	if (mpi_rank == mpi_root)
//		GGtsize = GGt.cols;
//
//
//	MPI_Bcast( & GGtsize, 1, esglobal_mpi, 0, environment->MPICommunicator);
//
//
//#if TIME_MEAS >= 1
//	double end = omp_get_wtime();
//	ESINFO(PROGRESS3) <<"CG Loop - Create GGt  - collect all matrices   - Runtime = " << ec1 - sc1 << " s";
//	ESINFO(PROGRESS3) <<"CG Loop - Create GGt  - GGt fact. processing   - Runtime = " << ep1 - sp1 << " s";
//	ESINFO(PROGRESS3) <<"CG Loop - Create GGt  - total = proc + comm    - Runtime = " << end - start << " s";
//#endif

}

void IterSolverBase::CreateGGt_Inv_old( SuperCluster & cluster )
{
	// To be removed
}

void IterSolverBase::CreateGGt_Inv( SuperCluster & cluster )
{

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Create_GGT_Inv, "Create_GGT_Inv", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	// temp variables
	vector < SparseMatrix > G_neighs   ( cluster.my_neighs.size() );
	vector < SparseMatrix > GGt_neighs ( cluster.my_neighs.size() );
	SparseMatrix 			G1t_l;
	SparseMatrix 			GGt_l;
	SparseMatrix 			GGt_Mat_tmp;
	SparseSolverMKL 		GGt_tmp;

    /* Numbers of processors, value of OMP_NUM_THREADS */
	int num_procs     = environment->PAR_NUM_THREADS;
	GGt_tmp.iparm[2]  = num_procs;

#ifdef READEX_LEVEL_2
	READEX_REGION_START(REG_Create_GGT_Inv_Exchange_local_GGt_MPI, "Create_GGT_Inv--Exchange_local_GGt_MPI", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 TimeEvent SaRGlocal("Exchange local G1 matrices to neighs. "); SaRGlocal.start();
	if (cluster.SYMMETRIC_SYSTEM)  {
		ExchangeMatrices(cluster.G1, G_neighs, cluster.my_neighs);
	} else {
		ExchangeMatrices(cluster.G2, G_neighs, cluster.my_neighs);
	}
	 SaRGlocal.end(); SaRGlocal.printStatMPI(); preproc_timing.addEvent(SaRGlocal);

#ifdef READEX_LEVEL_2
	READEX_REGION_STOP(REG_Create_GGT_Inv_Exchange_local_GGt_MPI);
#endif

	 TimeEvent Gt_l_trans("Local G1 matrix transpose to create Gt "); Gt_l_trans.start();
	if (cluster.USE_HFETI == 0) {
		cluster.G1.MatTranspose(G1t_l);
	}
	 Gt_l_trans.end(); Gt_l_trans.printStatMPI(); preproc_timing.addEvent(Gt_l_trans);

#ifdef READEX_LEVEL_2
	READEX_REGION_START(REG_Create_GGT_Inv_G1_MatMat_MatAdd, "Create_GGT_Inv--G1_MatMat+MatAdd", SCOREP_USER_REGION_TYPE_COMMON);
#endif

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

	int local_ker_size  = (int)cluster.G1.rows;
	int global_ker_size = 0;
	int global_GGt_size = 0;

	SEQ_VECTOR<int> global_ker_sizes;
	global_ker_sizes.resize(environment->MPIsize, 0);

	MPI_Exscan(&local_ker_size, &global_ker_size, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);
	MPI_Allgather(&global_ker_size, 1, MPI_INT, &global_ker_sizes[0],1, MPI_INT, environment->MPICommunicator);
	MPI_Allreduce(&local_ker_size, &global_GGt_size, 1, MPI_INT, MPI_SUM, environment->MPICommunicator);

	for (size_t i = 0; i < GGt_l.CSR_J_col_indices.size(); i++) {
		GGt_l.CSR_J_col_indices[i] += global_ker_size;
	}
	GGt_l.cols = global_GGt_size;


	 TimeEvent GGTNeighTime("G1t_local x G1_neigh MatMat(N-times) "); GGTNeighTime.start();
	#pragma omp parallel for
	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		GGt_neighs[neigh_i].MatMatT(G_neighs[neigh_i], cluster.G1);
		GGt_neighs[neigh_i].MatTranspose();
		eslocal inc = global_ker_sizes[cluster.my_neighs[neigh_i]];
		for (size_t i = 0; i < GGt_neighs[neigh_i].CSR_J_col_indices.size(); i++) {
			GGt_neighs[neigh_i].CSR_J_col_indices[i] += inc;
		}
		GGt_neighs[neigh_i].cols = global_GGt_size;
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

#ifdef READEX_LEVEL_2
	READEX_REGION_STOP(REG_Create_GGT_Inv_G1_MatMat_MatAdd);
#endif

#ifdef READEX_LEVEL_2
	READEX_REGION_START(REG_Create_GGT_Inv_CollectGGtPiecesToMaster, "Create_GGT_Inv--CollectGGtPiecesToMaster", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 // Collecting pieces of GGt from all clusters to master (MPI rank 0) node - using binary tree reduction
	 TimeEvent collectGGt_time("Collect GGt pieces to master"); 	collectGGt_time.start();
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

		MPI_Barrier(environment->MPICommunicator);

		GGt_l.Clear();
		count_cv_l += mpi_size/li;
		ESINFO(PROGRESS3) << "Collecting matrices G : " << count_cv_l <<" of " << mpi_size;
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

#ifdef READEX_LEVEL_2
	READEX_REGION_STOP(REG_Create_GGT_Inv_CollectGGtPiecesToMaster);
#endif

	// *** Calculating inverse GGt matrix in distributed fashion ***********************************************************
	// Create Sparse Direct solver for GGt
	if (mpi_rank == mpi_root) {
		GGt_tmp.msglvl = Info::report(LIBRARIES) ? 1 : 0;
	}

	 TimeEvent importGGt_time("Time to import GGt matrix into solver"); importGGt_time.start();
	GGt_Mat_tmp.mtype = cluster.mtype;
	GGt_tmp.ImportMatrix_wo_Copy (GGt_Mat_tmp);
	 importGGt_time.end(); importGGt_time.printStatMPI(); preproc_timing.addEvent(importGGt_time);

#ifdef READEX_LEVEL_2
	READEX_REGION_START(REG_Create_GGT_Inv_GGt_Factorization, "Create_GGT_Inv--GGt_Factorization", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 TimeEvent GGtFactor_time("GGT Factorization time"); GGtFactor_time.start();
	GGt_tmp.SetThreaded();
	std::stringstream ss;
	ss << "Create GGt_inv_dist-> rank: " << environment->MPIrank;
	GGt_tmp.Factorization(ss.str());
	 GGtFactor_time.end(); GGtFactor_time.printStatMPI(); preproc_timing.addEvent(GGtFactor_time);

#ifdef READEX_LEVEL_2
	READEX_REGION_STOP(REG_Create_GGT_Inv_GGt_Factorization);
#endif

	 TimeEvent GGT_rhs_time("Time to create InitialCondition for get GGTINV"); GGT_rhs_time.start();
	SEQ_VECTOR <double> rhs             (cluster.G1.rows * GGt_tmp.rows, 0.0);
	cluster.GGtinvM.dense_values.resize (cluster.G1.rows * GGt_tmp.rows, 0.0);
	for (eslocal i = 0; i < cluster.G1.rows; i++) {
		eslocal index = (GGt_tmp.rows * i) + i + global_ker_size;
		rhs[index] = 1;
	}
	 GGT_rhs_time.end(); GGT_rhs_time.printStatMPI(); preproc_timing.addEvent(GGT_rhs_time);

#ifdef READEX_LEVEL_2
	READEX_REGION_START(REG_Create_GGT_Inv_GGt_Solve, "Create_GGT_Inv--GGt_Solve", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 TimeEvent GGt_solve_time("Running solve to get stripe(s) of GGtINV"); GGt_solve_time.start();
	if (cluster.G1.rows > 0) {
		GGt_tmp.Solve(rhs, cluster.GGtinvM.dense_values, cluster.G1.rows);
	}

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

#ifdef READEX_LEVEL_2
	READEX_REGION_STOP(REG_Create_GGT_Inv_GGt_Solve);
#endif

	MKL_Set_Num_Threads(1);

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Create_GGT_Inv);
#endif
}


// *** Projector routines ************************************************
void IterSolverBase::Projector (TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
{

	ESINFO(GLOBAL_ERROR) << "Projector with factorized GGt matrix is not supported in current version";

	// TODO: Obsolete code - must be updated before used with current version
    // Code is commented

//	 time_eval.totalTime.start();
//	eslocal d_local_size = cluster.G1_comp.rows;
//	eslocal mpi_root     = 0;
//
//	SEQ_VECTOR<double> d_local( d_local_size );
//	SEQ_VECTOR<double> d_mpi  ( GGtsize );
//
//	 time_eval.timeEvents[0].start();
//	if ( output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 1)
//		d_local = x_in;
//	else
//		cluster.G1_comp.MatVec(x_in, d_local, 'N');
//	 time_eval.timeEvents[0].end();
//
//
//
//	 time_eval.timeEvents[1].start();
//	MPI_Gather(&d_local[0], d_local_size, MPI_DOUBLE,
//		&d_mpi[0], d_local_size, MPI_DOUBLE,
//		mpi_root, environment->MPICommunicator);
//	 time_eval.timeEvents[1].end();
//
//	time_eval.timeEvents[2].start();
//	if (mpi_rank == mpi_root ) {
//		GGt.Solve(d_mpi);				// t1 = Uc\(Lc\d);
//	}
//	time_eval.timeEvents[2].end();
//
//	time_eval.timeEvents[3].start();
//	MPI_Scatter( &d_mpi[0],      d_local_size, MPI_DOUBLE,
//		&d_local[0], d_local_size, MPI_DOUBLE,
//		mpi_root, environment->MPICommunicator);
//	time_eval.timeEvents[3].end();
//
//	if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 2) {
//		// for mu calculation
//		y_out = d_local;
//
//	} else {
//
//		time_eval.timeEvents[4].start();
//		//cluster.G1t_comp.MatVec(d_local, cluster.compressed_tmp, 'N'); // SUPER POZOR
//		cluster.G1_comp.MatVec(d_local, cluster.compressed_tmp, 'T');
//		time_eval.timeEvents[4].end();
//
//		time_eval.timeEvents[5].start();
//		All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
//		time_eval.timeEvents[5].end();
//
//		if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 0) {
//			#pragma omp parallel for
//			for (size_t i = 0; i < x_in.size(); i++)
//				y_out[i] = x_in[i] - y_out[i];
//		}
//
//	}
//
//	time_eval.totalTime.end();

}


/*
void IterSolverBase::Projector_l_inv_compG ( TimeEval & time_eval, Cluster & cluster, SparseMatrix & X_in, SparseMatrix & Y_out ) {

	Y_out.Clear();
	Y_out.type = 'G';
	Y_out.rows = X_in.rows;
	Y_out.cols = X_in.cols;
	Y_out.nnz  = X_in.dense_values.size();

	for (int i = 0; i < X_in.cols; i++) {

		std::vector<double> tmp_pr_in  (X_in.dense_values.begin() + i*X_in.rows,  X_in.dense_values.begin() + (i+1)*X_in.rows);
		std::vector<double> tmp_pr_out (tmp_pr_in.size(), 0 );

		Projector_l_inv_compG(timeEvalProj, cluster, tmp_pr_in, tmp_pr_out, 0  );

		Y_out.dense_values.insert(Y_out.dense_values.end(), tmp_pr_out.begin(), tmp_pr_out.end());
	}

}


void IterSolverBase::apply_A_l_Mat( TimeEval & time_eval, Cluster & cluster, SparseMatrix       & X_in, SparseMatrix       & Y_out) {


	Y_out.Clear();
	Y_out.type = 'G';
	Y_out.rows = X_in.rows;
	Y_out.cols = X_in.cols;
	Y_out.nnz  = X_in.dense_values.size();

	for (int i = 0; i < X_in.cols; i++) {

		std::vector<double> tmp_pr_in  (X_in.dense_values.begin() + i*X_in.rows,  X_in.dense_values.begin() + (i+1)*X_in.rows);
		std::vector<double> tmp_pr_out (tmp_pr_in.size(), 0 );

		apply_A_l_comp_dom_B(timeEvalAppa, cluster, tmp_pr_in, tmp_pr_out);

		Y_out.dense_values.insert(Y_out.dense_values.end(), tmp_pr_out.begin(), tmp_pr_out.end());
	}



}
*/

void IterSolverBase::Projector_Inv (TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
{

#ifdef READEX_LEVEL_1
	READEX_REGION_START(REG_Projector_Inv, "Projector_Inv", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	 time_eval.totalTime.start();
	eslocal d_local_size = cluster.G1_comp.rows;
	SEQ_VECTOR<double> d_local( d_local_size );
	SEQ_VECTOR<double> d_mpi  ( GGtsize );
	 time_eval.timeEvents[0].start();

	if (   output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 1
			||
		   output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 3)
	{
		d_local = x_in;
	} else {
		if (cluster.SYMMETRIC_SYSTEM) {
			if (cluster.G1_comp.nnz > 0) {
				cluster.G1_comp.MatVec(x_in, d_local, 'N');
			}
		} else {
			if (cluster.G2_comp.nnz > 0) {
				cluster.G2_comp.MatVec(x_in, d_local, 'N');
			}
		}
	}

	 time_eval.timeEvents[0].end();

#ifdef READEX_LEVEL_2
	READEX_REGION_START(REG_Projector_Inv_AllGather, "Projector_Inv--AllGather", SCOREP_USER_REGION_TYPE_COMMON);
#endif

	//TODO: Udelat poradne
	 time_eval.timeEvents[1].start();
	SEQ_VECTOR<int> ker_size_per_clusters(environment->MPIsize,0);
	MPI_Allgather(&d_local_size, 1, MPI_INT, &ker_size_per_clusters[0], 1, MPI_INT, environment->MPICommunicator );

	SEQ_VECTOR<int> displs (environment->MPIsize,0);
	displs[0] = 0;

	for (size_t i=1; i<displs.size(); ++i) {
		displs[i] = displs[i-1] + ker_size_per_clusters[i-1];
	}
	MPI_Allgatherv(&d_local[0], d_local_size, MPI_DOUBLE, &d_mpi[0], &ker_size_per_clusters[0], &displs[0], MPI_DOUBLE, environment->MPICommunicator);
	// TODO: END

	time_eval.timeEvents[1].end();

#ifdef READEX_LEVEL_2
	READEX_REGION_STOP(REG_Projector_Inv_AllGather);
#endif

#ifdef READEX_LEVEL_2
	READEX_REGION_START(REG_Projector_Inv_DenseMatVec, "Projector_Inv--DenseMatVec", SCOREP_USER_REGION_TYPE_COMMON);
#endif
	// TODO: END

	 time_eval.timeEvents[2].start();

	if (cluster.GGtinvM.cols != 0) {
		cluster.GGtinvM.DenseMatVec(d_mpi, d_local, 'T');
	}
	 time_eval.timeEvents[2].end();

#ifdef READEX_LEVEL_2
	READEX_REGION_STOP(REG_Projector_Inv_DenseMatVec);
#endif

	 time_eval.timeEvents[3].start();
	//MPI_Scatter( &d_mpi[0],      d_local_size, MPI_DOUBLE, &d_local[0], d_local_size, MPI_DOUBLE, mpi_root, environment->MPICommunicator);
	 time_eval.timeEvents[3].end();

	if (output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 2
		 ||
		output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0 == 3)
	{
		y_out = d_local; // for RBM amplitudes calculation
	} else {
		 time_eval.timeEvents[4].start();
		if (cluster.G1_comp.nnz > 0)
			cluster.G1_comp.MatVec(d_local, cluster.compressed_tmp, 'T');
		else
			std::fill (cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
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

#ifdef READEX_LEVEL_1
	READEX_REGION_STOP(REG_Projector_Inv);
#endif

}

void IterSolverBase::Projector_Inv_old (TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out, eslocal output_in_kerr_dim_2_input_in_kerr_dim_1_inputoutput_in_dual_dim_0) // eslocal mpi_rank, SparseSolverCPU & GGt,
{
	//TODO - To be completely removed
}

void IterSolverBase::Apply_Prec( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out )
{
	// implemeted in itersolvercpu, gpu and acc
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
		MPI_Send(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, environment->MPICommunicator);
		MPI_Send(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, environment->MPICommunicator );
		MPI_Send(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, environment->MPICommunicator );
		MPI_Send(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, environment->MPICommunicator );
#else
		MPI_Isend(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, environment->MPICommunicator, & request);
		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, environment->MPICommunicator, & request);
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, environment->MPICommunicator, & request);
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, environment->MPICommunicator, & request);
#endif

	}

	if (rank == dest_rank) {
		eslocal recv_par_buf[4];
		MPI_Recv(recv_par_buf, 4, esglobal_mpi, source_rank, param_tag, environment->MPICommunicator, & status);
		B_out.cols = recv_par_buf[0];
		B_out.rows = recv_par_buf[1];
		B_out.nnz  = recv_par_buf[2];
		B_out.type = recv_par_buf[3];

		B_out.CSR_I_row_indices.resize(B_out.rows + 1);
		B_out.CSR_J_col_indices.resize(B_out.nnz);
		B_out.CSR_V_values.     resize(B_out.nnz);

		MPI_Recv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, esglobal_mpi,    source_rank, I_row_tag, environment->MPICommunicator, & status );
		MPI_Recv(&B_out.CSR_J_col_indices[0], B_out.nnz,      esglobal_mpi,    source_rank, J_col_tag, environment->MPICommunicator, & status );
		MPI_Recv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE, source_rank, V_val_tag, environment->MPICommunicator, & status );
	}

#ifdef WIN32
	MPI_Barrier(environment->MPICommunicator);
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

		MPI_Isend(send_par_buf, 		   				  4, 	esglobal_mpi, 	dest_rank, tag, environment->MPICommunicator, &request[0] );
		if (A_in.nnz > 0) {
			MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, 	esglobal_mpi, 	dest_rank, tag, environment->MPICommunicator, &request[1] );
			MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      	esglobal_mpi, 	dest_rank, tag, environment->MPICommunicator, &request[2] );
			MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   		MPI_DOUBLE, 	dest_rank, tag, environment->MPICommunicator, &request[3] );
			MPI_Waitall( 4 , &request[0], MPI_STATUSES_IGNORE);
		} else {
			// Empty matrix
			MPI_Waitall( 1 , &request[0], MPI_STATUSES_IGNORE);
		}
	}

	if (rank == dest_rank) {

		MPI_Status status;
		SEQ_VECTOR < MPI_Request > request ( 3 );
		eslocal recv_par_buf[4];

		MPI_Recv(recv_par_buf, 4, esglobal_mpi, source_rank, tag, environment->MPICommunicator, & status);
		B_out.cols = recv_par_buf[0];
		B_out.rows = recv_par_buf[1];
		B_out.nnz  = recv_par_buf[2];
		B_out.type = recv_par_buf[3];

		if (B_out.nnz > 0) {
			B_out.CSR_I_row_indices.resize(B_out.rows + 1);
			B_out.CSR_J_col_indices.resize(B_out.nnz);
			B_out.CSR_V_values.     resize(B_out.nnz);

			MPI_Irecv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, esglobal_mpi,    source_rank, tag, environment->MPICommunicator, &request[0] );
			MPI_Irecv(&B_out.CSR_J_col_indices[0], B_out.nnz,      esglobal_mpi,    source_rank, tag, environment->MPICommunicator, &request[1] );
			MPI_Irecv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE,      source_rank, tag, environment->MPICommunicator, &request[2] );

			MPI_Waitall( 3 , &request[0], MPI_STATUSES_IGNORE);
		}
	}

}


void   SendMatrix ( SparseMatrix & A_in, eslocal dest_rank ) {

	int rank;
	MPI_Comm_rank (environment->MPICommunicator, &rank);


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
//		MPI_Send(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, environment->MPICommunicator);
//		MPI_Send(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, environment->MPICommunicator );
//		MPI_Send(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, environment->MPICommunicator );
//		MPI_Send(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, environment->MPICommunicator );
//#else
		MPI_Isend(send_par_buf, 4, esglobal_mpi, dest_rank, param_tag, environment->MPICommunicator, & request);
		MPI_Isend(&A_in.CSR_I_row_indices[0], A_in.rows + 1, esglobal_mpi, dest_rank, I_row_tag, environment->MPICommunicator, & request);
		MPI_Isend(&A_in.CSR_J_col_indices[0], A_in.nnz,      esglobal_mpi, dest_rank, J_col_tag, environment->MPICommunicator, & request);
		MPI_Isend(&A_in.CSR_V_values[0],      A_in.nnz,   MPI_DOUBLE, dest_rank, V_val_tag, environment->MPICommunicator, & request);
//#endif

#ifdef WIN32
//	MPI_Barrier(environment->MPICommunicator);
#endif
}

void   RecvMatrix ( SparseMatrix & B_out, eslocal source_rank) {

	int rank;
	MPI_Comm_rank (environment->MPICommunicator, &rank);


	eslocal param_tag = 1;
	eslocal I_row_tag = 2;
	eslocal J_col_tag = 3;
	eslocal V_val_tag = 4;

	MPI_Status status;

	eslocal recv_par_buf[4];
	MPI_Recv(recv_par_buf, 4, esglobal_mpi, source_rank, param_tag, environment->MPICommunicator, & status);
	B_out.cols = recv_par_buf[0];
	B_out.rows = recv_par_buf[1];
	B_out.nnz  = recv_par_buf[2];
	B_out.type = recv_par_buf[3];

	B_out.CSR_I_row_indices.resize(B_out.rows + 1);
	B_out.CSR_J_col_indices.resize(B_out.nnz);
	B_out.CSR_V_values.     resize(B_out.nnz);

	MPI_Recv(&B_out.CSR_I_row_indices[0], B_out.rows + 1, esglobal_mpi,    source_rank, I_row_tag, environment->MPICommunicator, & status );
	MPI_Recv(&B_out.CSR_J_col_indices[0], B_out.nnz,      esglobal_mpi,    source_rank, J_col_tag, environment->MPICommunicator, & status );
	MPI_Recv(&B_out.CSR_V_values[0],      B_out.nnz,      MPI_DOUBLE, source_rank, V_val_tag, environment->MPICommunicator, & status );


#ifdef WIN32
	//MPI_Barrier(environment->MPICommunicator);
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

		MPI_Isend(&send_par_buf[neigh_i][0],  4, 				esglobal_mpi, 	dest_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 0] );

	}


	for (size_t neigh_i = 0; neigh_i < neighbor_ranks.size(); neigh_i++ )
	{
		recv_par_buf[neigh_i].resize(4);
		eslocal source_rank = neighbor_ranks[neigh_i];

		MPI_Recv(&recv_par_buf[neigh_i][0], 4, esglobal_mpi, source_rank, tag, environment->MPICommunicator, MPI_STATUS_IGNORE);

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

		if (A_in.rows) {
			MPI_Isend(A_in.CSR_I_row_indices.data(), A_in.rows + 1, 	esglobal_mpi, 	dest_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 1] );
		} else {
			MPI_Isend(A_in.CSR_I_row_indices.data(),             0, 	esglobal_mpi, 	dest_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 1] );
		}
		MPI_Isend(A_in.CSR_J_col_indices.data(), A_in.nnz,      	esglobal_mpi, 	dest_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 2] );
		MPI_Isend(A_in.CSR_V_values.data(),      A_in.nnz,   		MPI_DOUBLE, 	dest_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 3] );

	}


	for (size_t neigh_i = 0; neigh_i < neighbor_ranks.size(); neigh_i++ )
	{
		eslocal source_rank = neighbor_ranks[neigh_i];
		eslocal tag = 1;

		if (B_out[neigh_i].rows) {
			MPI_Irecv(&B_out[neigh_i].CSR_I_row_indices[0], B_out[neigh_i].rows + 1, esglobal_mpi,    source_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 4] );
		} else {
			MPI_Irecv(&B_out[neigh_i].CSR_I_row_indices[0],                       0, esglobal_mpi,    source_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 4] );
		}
		MPI_Irecv(&B_out[neigh_i].CSR_J_col_indices[0], B_out[neigh_i].nnz,      esglobal_mpi,    source_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 5] );
		MPI_Irecv(&B_out[neigh_i].CSR_V_values[0],      B_out[neigh_i].nnz,      MPI_DOUBLE,      source_rank, tag, environment->MPICommunicator, &request[7 * neigh_i + 6] );
	}

	MPI_Waitall(7 * neighbor_ranks.size(), &request[0], MPI_STATUSES_IGNORE);

}

void   BcastMatrix ( eslocal rank, eslocal mpi_root, eslocal source_rank, SparseMatrix & A) {

	eslocal send_par_buf[4];

	if (rank == source_rank) {
		send_par_buf[0] = A.cols; send_par_buf[1] = A.rows; send_par_buf[2] = A.nnz; send_par_buf[3] = A.type;
	}

	MPI_Bcast(send_par_buf, 4, esglobal_mpi, source_rank, environment->MPICommunicator);

	if (rank != source_rank) {
		A.cols = send_par_buf[0]; A.rows = send_par_buf[1]; A.nnz  = send_par_buf[2]; A.type = send_par_buf[3];
		A.CSR_I_row_indices.resize(A.rows + 1);
		A.CSR_J_col_indices.resize(A.nnz);
		A.CSR_V_values.     resize(A.nnz);
	}

	if (A.nnz > 0) {

		MPI_Bcast(&A.CSR_I_row_indices[0], A.rows + 1, esglobal_mpi, source_rank, environment->MPICommunicator);
		MPI_Bcast(&A.CSR_J_col_indices[0], A.nnz,      esglobal_mpi, source_rank, environment->MPICommunicator);
		MPI_Bcast(&A.CSR_V_values[0],      A.nnz,   MPI_DOUBLE, source_rank, environment->MPICommunicator);

	}
}

//void   All_Reduce_lambdas_compB2( Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out )
//{
//
//	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
//		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
//			cluster.my_comm_lambdas[i][j] = x_in[cluster.my_comm_lambdas_indices_comp[i][j]];
//		}
//	}
//
//
//	MPI_Request * mpi_req  = new MPI_Request [cluster.my_neighs.size()];
//	MPI_Status  * mpi_stat = new MPI_Status  [cluster.my_neighs.size()];
//
//	cluster.iter_cnt_comm++;
//	eslocal tag = cluster.iter_cnt_comm;
//
//	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
//		MPI_Sendrecv(
//			&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,
//			&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,
//			environment->MPICommunicator, &mpi_stat[neigh_i] );
//	}
//
//	//for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
//	//	size_t b_size = cluster.my_comm_lambdas[neigh_i].size();
//	//	MPI_Isend(&b_size,                              1                                      , esglobal_mpi   , cluster.my_neighs[neigh_i], tag + 100, environment->MPICommunicator, &mpi_req[neigh_i] );
//	//	MPI_Isend(&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag,       environment->MPICommunicator, &mpi_req[neigh_i] );
//	//
//	//}
//
//	//for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
//	//	size_t r_size = 0;
//	//	MPI_Recv(&r_size                             ,                                       1, esglobal_mpi   , cluster.my_neighs[neigh_i], tag + 100, environment->MPICommunicator, &mpi_stat[neigh_i] );
//	//	if (r_size != cluster.my_recv_lambdas[neigh_i].size()) cout << "Error - different buffer size " << endl;
//	//	MPI_Recv(&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag      , environment->MPICommunicator, &mpi_stat[neigh_i] );
//	//}
//
//#ifdef XE6
//	MPI_Barrier(environment->MPICommunicator);
//#endif
//
//#ifdef WIN32
//	MPI_Barrier(environment->MPICommunicator);
//#endif
//
//	delete [] mpi_req;
//	delete [] mpi_stat;
//
//	y_out = x_in; // POZOR pozor
//	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
//		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
//			y_out[cluster.my_comm_lambdas_indices_comp[i][j]] += cluster.my_recv_lambdas[i][j];
//		}
//	}
//
//
//
//}


void   All_Reduce_lambdas_compB( SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out )
{

	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			cluster.my_comm_lambdas[i][j] = x_in[cluster.my_comm_lambdas_indices_comp[i][j]];
		}
	}

	SEQ_VECTOR < MPI_Request > request ( 2 * cluster.my_neighs.size() );

	//cluster.iter_cnt_comm++;
	eslocal tag = 1;

	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		MPI_Isend(
			&cluster.my_comm_lambdas[neigh_i][0], cluster.my_comm_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag, environment->MPICommunicator, &request[ 0                        + neigh_i] );
	}

	for (size_t neigh_i = 0; neigh_i < cluster.my_neighs.size(); neigh_i++ ) {
		MPI_Irecv(
			&cluster.my_recv_lambdas[neigh_i][0], cluster.my_recv_lambdas[neigh_i].size(), MPI_DOUBLE, cluster.my_neighs[neigh_i], tag, environment->MPICommunicator, &request[ cluster.my_neighs.size() + neigh_i] );
	}

	MPI_Waitall( 2 * cluster.my_neighs.size(), &request[0], MPI_STATUSES_IGNORE);

	y_out = x_in; // TODO: POZOR pozor
	for (size_t i = 0; i < cluster.my_comm_lambdas_indices_comp.size(); i++) {
		for (size_t j = 0; j < cluster.my_comm_lambdas_indices_comp[i].size(); j++) {
			y_out[cluster.my_comm_lambdas_indices_comp[i][j]] += cluster.my_recv_lambdas[i][j];
		}
	}

}

void   compress_lambda_vector  ( SuperCluster & cluster, SEQ_VECTOR <double> & decompressed_vec_lambda)
{
	//compress vector for CG in main loop
	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[i] = decompressed_vec_lambda[cluster.my_lamdas_indices[i]];

	decompressed_vec_lambda.resize(cluster.my_lamdas_indices.size());
}

void   decompress_lambda_vector( SuperCluster & cluster, SEQ_VECTOR <double> & compressed_vec_lambda)
{
	SEQ_VECTOR <double> decompressed_vec_lambda (cluster.domains[0]->B1.rows,0); //TODO : needs fix

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)
		decompressed_vec_lambda[cluster.my_lamdas_indices[i]] = compressed_vec_lambda[i];

	compressed_vec_lambda = decompressed_vec_lambda;
}

double parallel_norm_compressed( SuperCluster & cluster, SEQ_VECTOR<double> & input_vector )
{

	double wl = 0; double wg = 0;

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)
		wl = wl + (input_vector[i] * input_vector[i] * cluster.my_lamdas_ddot_filter[i]);

	MPI_Allreduce( &wl, &wg, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
	double norm_l = sqrt(wg);

	return norm_l;
}


double parallel_ddot_compressed_double( SuperCluster & cluster, double * input_vector1, double * input_vector2 )
{
	double a1 = 0; double a1g = 0;

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
		a1 = a1 + (input_vector1[i] * input_vector2[i] * cluster.my_lamdas_ddot_filter[i]);
	}

	MPI_Allreduce( &a1, &a1g, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);

	return a1g;
}


double parallel_ddot_compressed( SuperCluster & cluster, SEQ_VECTOR<double> & input_vector1, SEQ_VECTOR<double> & input_vector2 )
{
	double a1 = 0; double a1g = 0;

	for (size_t i = 0; i < cluster.my_lamdas_indices.size(); i++)  {
		a1 = a1 + (input_vector1[i] * input_vector2[i] * cluster.my_lamdas_ddot_filter[i]);
	}

	MPI_Allreduce( &a1, &a1g, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);

	return a1g;
}

void   parallel_ddot_compressed_non_blocking( SuperCluster & cluster,
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
	MPI_Barrier(environment->MPICommunicator);
	MPI_Allreduce( &send_buf[0], &output[0], 3, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
#else
#ifdef USE_MPI_3
	MPI_Iallreduce( &send_buf[0], &output[0], 3, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator, mpi_req);
#else
	MPI_Allreduce( &send_buf[0], &output[0], 3, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
#endif
#endif

}


void   parallel_ddot_compressed_non_blocking( SuperCluster & cluster,
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
	MPI_Barrier(environment->MPICommunicator);
	MPI_Allreduce( &send_buf[0], &output[0], 2, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
#else
#ifdef USE_MPI_3
	MPI_Iallreduce( &send_buf[0], &output[0], 2, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator, mpi_req);
#else
	MPI_Allreduce( &send_buf[0], &output[0], 2, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
#endif
#endif

}

}

// *** END - Communication layer ****************************************
// **********************************************************************





