#include "itersolverGPU.h"

using namespace espreso;

// *** Action of K+ routines *********************************************

void IterSolverGPU::apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {
       time_eval.totalTime.start();

	if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {
		time_eval.timeEvents[0].start();

		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {

			// *** Part 1 - prepare vectors for FETI operator with SC
			//     cluster.domains[d].compressed_tmp2 for CPU
			//     cluster.domains[d].cuda_pinned_buff - for GPU - double precision
			//     cluster.domains[d].cuda_pinned_buff_fl[i] - for GPU - single precision

			if (cluster.domains[d].isOnACC == 1) {
				if (!cluster.domains[d].B1Kplus.USE_FLOAT) {
					for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.domains[d].cuda_pinned_buff[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
						cluster.domains[d].compressed_tmp2[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
					}
				} else {
					for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.domains[d].cuda_pinned_buff_fl[i] = (float) x_in[ cluster.domains[d].lambda_map_sub_local[i]];
						cluster.domains[d].compressed_tmp2[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
					}
				}
			} else {
				for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
					cluster.domains[d].compressed_tmp2[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
				}
			}

			// *** Part 2 - Prepare data for HTFETI operator (multKplusGLobal_Kinv
			//     on CPU - cluster.x_prim_cluster1[d]
			cluster.domains[d].B1_comp_dom.MatVec (cluster.domains[d].compressed_tmp2, cluster.x_prim_cluster1[d], 'T');
		}

//		// *** Part 3 - execute FETI SC operator - B1Kplus.DenseMatVec
//		// Note: DenseMatVecCUDA_wo_Copy_start is non-blocking operation - just schedule transfers and execution of kernels and goes forward
//		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
//		  	if (cluster.domains[d].isOnACC == 1) {
//				if (!cluster.domains[d].B1Kplus.USE_FLOAT) {
//					//cilk_spawn
//					cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start( cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff,'N',0 );
//				} else {
//					//cilk_spawn
//					cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start_fl( cluster.domains[d].cuda_pinned_buff_fl, cluster.domains[d].cuda_pinned_buff_fl,'N',0 );
//				}
//		  	}
//		//}
//
//		//cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
//		  	if (cluster.domains[d].isOnACC == 0) {
//				// Automatic fall-back to CPU for sub-domains, which did not fit GPU memory
//				// Note: DenseMatVec - this is a blocking operation - it waits till its finished
//		  		//cilk_spawn
//				cluster.domains[d].B1Kplus.DenseMatVec (cluster.domains[d].compressed_tmp2, cluster.domains[d].compressed_tmp);
//			}
//		}

		cilk_spawn apply_A_feti_SC   ( cluster );

		time_eval.timeEvents[0].end();

		// *** Part 4 - Execute HTFETI operator
		time_eval.timeEvents[1].start();
		//cilk_spawn
		apply_A_htfeti_SC ( cluster );
		cilk_sync;
		time_eval.timeEvents[1].end();

		time_eval.timeEvents[2].start();

		// *** Part 5 - Finalize transfers of the result of the FETI SC operator from GPU back to CPU
		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d].isOnACC == 1) {
				cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_sync ( );
			}
		}

		// Part 6 - Lambda values, results of the FETI SC operator, are combined back into single lambda vector per cluster - cluster.compressed_tmp
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d].isOnACC == 1) {
				if (!cluster.domains[d].B1Kplus.USE_FLOAT) {
					for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].cuda_pinned_buff[i];
					}
				} else {
					for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += (double)cluster.domains[d].cuda_pinned_buff_fl[i];
					}
				}
			} else {
					for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
					}
			}
		}

		// *** Part 7 - calculate lambda values from the results of the HTFETI operator (results are in cluster.x_prim_cluster1)
		//              using cluster.domains[d].B1_comp_dom gluing matrix and store them in the temp. vector y_out_tmp.
		//				Final result is combination of lambdas achieved by the FETI operator (cluster.compressed_tmp) and HTFETI operators (y_out_tmp)
		SEQ_VECTOR < double > y_out_tmp;
		for (eslocal d = 0; d < cluster.domains.size(); d++) {
			// Calculate lambdas from the HTFETI operator
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0);

			// Sums lambdas from FETI and HFETI operator
			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		time_eval.timeEvents[2].end();

	}
	time_eval.totalTime.start();




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

//#ifdef CUDA
			if (cluster.domains[d].isOnACC == 1) {
				cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy(x_in_tmp, cluster.domains[d].compressed_tmp,'N',0);
			} else {
				cluster.domains[d].B1Kplus.DenseMatVec            (x_in_tmp, cluster.domains[d].compressed_tmp);
			}
//#else
//			cluster.domains[d].B1Kplus.DenseMatVec            (x_in_tmp, cluster.domains[d].compressed_tmp);
//#endif
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







	if (cluster.USE_KINV == 0) {

		 time_eval.timeEvents[0].start();
		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
			for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
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


	// *** Combine update lambdas among neighbors - shared by all methos above

	time_eval.timeEvents[3].start();
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[3].end();

	time_eval.totalTime.end();

}


	void IterSolverGPU::apply_A_feti_SC   ( Cluster & cluster )
	{

		cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
		  	if (cluster.domains[d].isOnACC == 1) {
				if (!cluster.domains[d].B1Kplus.USE_FLOAT) {
					//cilk_spawn
					cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start( cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff,'N',0 );
				} else {
					//cilk_spawn
					cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start_fl( cluster.domains[d].cuda_pinned_buff_fl, cluster.domains[d].cuda_pinned_buff_fl,'N',0 );
				}
		  	}
		//}

		//cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
		  	if (!config::solver::COMBINE_SC_AND_SPDS) {
				if (cluster.domains[d].isOnACC == 0) {
					// Automatic fall-back to CPU for sub-domains Schur complements, which did not fit GPU memory
					// Note: DenseMatVec - this is a blocking operation - it waits till its finished
					//cilk_spawn
					cluster.domains[d].B1Kplus.DenseMatVec (cluster.domains[d].compressed_tmp2, cluster.domains[d].compressed_tmp);
				}
		  	}
		}
	}


	void IterSolverGPU::apply_A_htfeti_SC ( Cluster & cluster )
	{
		//cluster.multKplusGlobal_Kinv( cluster.x_prim_cluster1 );
		cluster.multKplusGlobal_GPU( cluster.x_prim_cluster1 );
	}
