#include "itersolverGPU.h"
#include "wrappers/nvtx/w.nvtx.h"

#include <iostream>

using namespace espreso;

//#define SHARE_SC

// *** Action of K+ routines *********************************************

void IterSolverGPU::apply_A_l_comp_dom_B( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {
       time_eval.totalTime.start();

       //ESINFO(ERROR) << "Implement apply_A_l_comp_dom_B.";
       // TODO: implement
       /*
	if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {
		time_eval.timeEvents[0].start();

		#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++) {

			// *** Part 1 - prepare vectors for FETI operator with SC
			//     cluster.domains[d].compressed_tmp2 for CPU
			//     cluster.domains[d].cuda_pinned_buff - for GPU - double precision
			//     cluster.domains[d].cuda_pinned_buff_fl[i] - for GPU - single precision

			if (cluster.domains[d].is_on_acc == 1) {
				if (!cluster.domains[d].B1Kplus.USE_FLOAT) {
					for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.domains[d].cuda_pinned_buff[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
						cluster.domains[d].compressed_tmp2[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
					}
				} else {
					for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.domains[d].cuda_pinned_buff_fl[i] = (float) x_in[ cluster.domains[d].lambda_map_sub_local[i]];
						cluster.domains[d].compressed_tmp2[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
					}
				}
			} else {
				for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
					cluster.domains[d].compressed_tmp2[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
				}
			}

			// *** Part 2 - Prepare data for HTFETI operator (multKplusGLobal_Kinv
			//     on CPU - cluster.x_prim_cluster1[d]
			cluster.domains[d].B1_comp_dom.MatVec (cluster.domains[d].compressed_tmp2, cluster.x_prim_cluster1[d], 'T');
		}

//		// *** Part 3 - execute FETI SC operator - B1Kplus.DenseMatVec
//		// Note: DenseMatVecCUDA_wo_Copy_start is non-blocking operation - just schedule transfers and execution of kernels and goes forward
//		cilk_for (esint d = 0; d < cluster.domains.size(); d++) {
//		  	if (cluster.domains[d].is_on_acc == 1) {
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
//		//cilk_for (esint d = 0; d < cluster.domains.size(); d++) {
//		  	if (cluster.domains[d].is_on_acc == 0) {
//				// Automatic fall-back to CPU for sub-domains, which did not fit GPU memory
//				// Note: DenseMatVec - this is a blocking operation - it waits till its finished
//		  		//cilk_spawn
//				cluster.domains[d].B1Kplus.DenseMatVec (cluster.domains[d].compressed_tmp2, cluster.domains[d].compressed_tmp);
//			}
//		}

		cilk_spawn apply_A_FETI_SC_forHFETI   ( cluster );

		time_eval.timeEvents[0].end();

		// *** Part 4 - Execute HTFETI operator
		time_eval.timeEvents[1].start();
		//cilk_spawn
		apply_A_htfeti_SC ( cluster );
		cilk_sync;
		time_eval.timeEvents[1].end();

		time_eval.timeEvents[2].start();

		// *** Part 5 - Finalize transfers of the result of the FETI SC operator from GPU back to CPU
		#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d].is_on_acc == 1) {
				cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_sync ( );
			}
		}

		// Part 6 - Lambda values, results of the FETI SC operator, are combined back into single lambda vector per cluster - cluster.compressed_tmp
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		for (esint d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d].is_on_acc == 1) {
				if (!cluster.domains[d].B1Kplus.USE_FLOAT) {
					for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].cuda_pinned_buff[i];
					}
				} else {
					for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += (double)cluster.domains[d].cuda_pinned_buff_fl[i];
					}
				}
			} else {
					for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
					}
			}
		}

		// *** Part 7 - calculate lambda values from the results of the HTFETI operator (results are in cluster.x_prim_cluster1)
		//              using cluster.domains[d].B1_comp_dom gluing matrix and store them in the temp. vector y_out_tmp.
		//				Final result is combination of lambdas achieved by the FETI operator (cluster.compressed_tmp) and HTFETI operators (y_out_tmp)
		SEQ_VECTOR < double > y_out_tmp;
		for (esint d = 0; d < cluster.domains.size(); d++) {
			// Calculate lambdas from the HTFETI operator
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0);

			// Sums lambdas from FETI and HFETI operator
			for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		time_eval.timeEvents[2].end();

	}
*/
	if (cluster.USE_HFETI == 0) {	//if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 0) {	
		
		// time_eval.timeEvents[0].start();
		//// POZOR - jen pro porovnani vykonu CPU a GPU
		//cilk_for (esint d = 0; d < cluster.domains.size(); d++) {
		//	SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
		//	for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
		//		x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
		//	cluster.domains[d].B1Kplus.DenseMatVec(x_in_tmp, cluster.domains[d].compressed_tmp);
		//}
		//// POZOR - jen pro porovnani vykonu CPU a GPU
		// time_eval.timeEvents[0].end();

		// apply_A_FETI_SC_forFETI (cluster, x_in);
	
		// Vsechny domeny ktere se vlezly na GPU se zpracuji zde 
		// 1. cast - vyzobani dat do pinovanych bufferu 
		time_eval.timeEvents[0].start();
		#pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
				for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
					cluster.domains[d]->cuda_pinned_buff[i] = x_in[ cluster.domains[d]->lambda_map_sub_local[i]];
				}
			}
		}
		time_eval.timeEvents[0].end();

		// 2. cast - nastartovani asynchroniho kopirovani vetktoru na GPU, DGEMV na GPU a kopirovani vektoru zpet na CPU 
		time_eval.timeEvents[1].start();
		PUSH_RANGE("apply_A_l_comp_dom_B GEMV", 0)
		// #pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
        		cluster.domains[d]->B1Kplus.DenseMatVecCUDA_wo_Copy_start(cluster.domains[d]->cuda_pinned_buff.data(), cluster.domains[d]->cuda_pinned_buff.data(),'N',0);
			}
		}

		// 2. cast - vsechny domeny, ktere se nevlezly na GPU se zpracuji na CPU - provede se zatimco se na GPU asynchrone provadi zpracovani domen lezicich na GPU 
		#pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d]->B1Kplus.is_on_acc == 0) {
				
				SEQ_VECTOR < double > x_in_tmp;
				x_in_tmp.resize( cluster.domains[d]->B1_comp_dom.rows );
				
				for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
					x_in_tmp[i] = x_in[ cluster.domains[d]->lambda_map_sub_local[i]];
				}

				if (cluster.USE_KINV == 0 || configuration.combine_sc_and_spds ) {
					cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'T');
					cluster.domains[d]->multKplusLocal(*cluster.x_prim_cluster1[d]);
					cluster.domains[d]->B1_comp_dom.MatVec (*cluster.x_prim_cluster1[d], cluster.domains[d]->compressed_tmp, 'N', 0, 0, 0.0);
				} else {
					cluster.domains[d]->B1Kplus.DenseMatVec (x_in_tmp, cluster.domains[d]->compressed_tmp);
				}

				// *** Original version 
				//
				// if (!configuration.combine_sc_and_spds) {
				// 	cluster.domains[d]->B1Kplus.DenseMatVec (x_in_tmp, cluster.domains[d]->compressed_tmp);
				// } else {
				// 	cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'T');
				// 	cluster.domains[d]->multKplusLocal(*cluster.x_prim_cluster1[d]);
				// 	cluster.domains[d]->B1_comp_dom.MatVec (*cluster.x_prim_cluster1[d], cluster.domains[d]->compressed_tmp, 'N', 0, 0, 0.0);
				// }
			}
		}

		// 3. cast - synchrizace exekuce na GPU
		#pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
				cluster.domains[d]->B1Kplus.DenseMatVecCUDA_wo_Copy_sync ( );
			}
		}
		POP_RANGE //apply_A_l_comp_dom_B GEMV

		time_eval.timeEvents[1].end();


		time_eval.timeEvents[2].start();

		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);

		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
				cluster.domains[d]->B1Kplus.DenseMatVecCUDA_wo_Copy_sync ( );
				for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
					cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += cluster.domains[d]->cuda_pinned_buff[i];
				}
			} else {
				for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++)
				{
					cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += cluster.domains[d]->compressed_tmp[i];
				}
			}
		}
		time_eval.timeEvents[2].end();
	}

/*

	if (cluster.USE_KINV == 0 && cluster.USE_HFETI == 1) {

		 time_eval.timeEvents[0].start();
		#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
			for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
		}
		 time_eval.timeEvents[0].end();

		 time_eval.timeEvents[1].start();

		int method = 0;

		if (method == 0) {
			cluster.multKplusGlobal_l(cluster.x_prim_cluster1);
		}

		if (method == 1) {
				cluster.multKplus_HF(cluster.x_prim_cluster1);

				#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++) {
					cluster.domains[d].multKplusLocal( cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d] );
				}
		}

		if (method == 2) {
			cilk_spawn cluster.multKplus_HF_SC (cluster.x_prim_cluster1, cluster.x_prim_cluster3);

			#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++) {
				cluster.domains[d].multKplusLocal( cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d] );
			}

			cilk_sync;

			#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++) {
				for (int i = 0; i < cluster.domains[d].K.rows; i++) {
					cluster.x_prim_cluster1[d][i] = cluster.x_prim_cluster2[d][i] + cluster.x_prim_cluster3[d][i];
				}
			}


		}

		 time_eval.timeEvents[1].end();


		 time_eval.timeEvents[2].start();
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (esint d = 0; d < cluster.domains.size(); d++) {
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		 time_eval.timeEvents[2].end();


	}


	if (cluster.USE_KINV == 0 && cluster.USE_HFETI == 0) {

		 time_eval.timeEvents[0].start();
		#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
			for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
		}
		 time_eval.timeEvents[0].end();

		 time_eval.timeEvents[1].start();
		if (cluster.USE_HFETI == 0) {
			#pragma omp parallel for
for (esint d = 0; d < cluster.domains.size(); d++)
				cluster.domains[d].multKplusLocal(cluster.x_prim_cluster1[d]);
		} else {
			cluster.multKplusGlobal_l(cluster.x_prim_cluster1);
		}
		 time_eval.timeEvents[1].end();


		 time_eval.timeEvents[2].start();
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
		SEQ_VECTOR < double > y_out_tmp;
		for (esint d = 0; d < cluster.domains.size(); d++) {
			y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

			for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
				cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
		}
		 time_eval.timeEvents[2].end();

	}
*/

	// *** Combine update lambdas among neighbors - shared by all methos above

	time_eval.timeEvents[3].start();
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	time_eval.timeEvents[3].end();

	time_eval.totalTime.end();

}



void IterSolverGPU::apply_A_FETI_SC_forFETI   ( SuperCluster & cluster, SEQ_VECTOR<double> & x_in ) {
		//cilk_
		// TODO_GPU - spustit v OpeMP smycce - mozna ne nutne "DenseMatVecCUDA_wo_Copy_start" ale urcite tu vyzobavaci smycku 
	#pragma omp parallel for
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
			for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
				cluster.domains[d]->cuda_pinned_buff[i] = x_in[ cluster.domains[d]->lambda_map_sub_local[i]];
			}
		}
	}
	// PUSH_RANGE("apply_A_FETI_SC_forFETI GEMV", 1) - multiplication sync. not found (POP_RANGE location)
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
#ifdef SHARE_SC
				cluster.domains[d]->B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d]->cuda_pinned_buff, cluster.domains[d]->cuda_pinned_buff,'N',0);
#else
				// TODO_GPU - tato funkce je ve SparseMatrix.cpp 
                                cluster.domains[d]->B1Kplus.DenseMatVecCUDA_wo_Copy_start(cluster.domains[d]->cuda_pinned_buff.data(), cluster.domains[d]->cuda_pinned_buff.data(),'N',0);
#endif
			}


//			if (cluster.domains[d].is_on_acc == 1 && d%2 == 0) {
//				for (esint i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
//					cluster.domains[d].cuda_pinned_buff[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
//				}
//
////				cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start(cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff,'N',0);
//
//				esint sc1_rows = cluster.domains[d].B1Kplus.rows;
//				esint sc2_rows = 0;
//
//				if(d+1 < cluster.domains.size()) {
//					for (esint i = 0; i < cluster.domains[d+1].lambda_map_sub_local.size(); i++) {
//						cluster.domains[d+1].cuda_pinned_buff[i] = x_in[ cluster.domains[d+1].lambda_map_sub_local[i]];
//					}
//					sc2_rows = cluster.domains[d+1].B1Kplus.rows;
//				}
//
//				if(sc1_rows > sc2_rows) {
//					cluster.domains[d].B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff, 'N', 0, sc1_rows, sc1_rows, 0, 'U');
//					if(sc2_rows > 0) {
//						cluster.domains[d].B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d+1].cuda_pinned_buff, cluster.domains[d+1].cuda_pinned_buff, 'N', 0, sc2_rows, sc1_rows, 1, 'L');
//					}
//				} else if(sc1_rows == sc2_rows) {
//					cluster.domains[d].B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff, 'N', 0, sc1_rows, sc1_rows, sc1_rows, 'U');
//					cluster.domains[d].B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d+1].cuda_pinned_buff, cluster.domains[d+1].cuda_pinned_buff, 'N', 0, sc2_rows, sc2_rows, 0, 'L');
//				} else { // sc1_rows < sc2_rows
//					cluster.domains[d].B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d].cuda_pinned_buff, cluster.domains[d].cuda_pinned_buff, 'N', 0, sc1_rows, sc2_rows, 1, 'L');
//					cluster.domains[d].B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d+1].cuda_pinned_buff, cluster.domains[d+1].cuda_pinned_buff, 'N', 0, sc2_rows, sc2_rows, 0, 'U');
//				}
//			}
		}
	}

void IterSolverGPU::apply_A_FETI_SC_forHFETI   ( Cluster & cluster )
	{
		//cilk_
		// PUSH_RANGE("apply_A_FETI_SC_forHFETI GEMV", 2)  - multiplication sync. not found (POP_RANGE location)
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d].B1Kplus.is_on_acc == 1) {
				if (!cluster.domains[d].B1Kplus.USE_FLOAT) {
					//cilk_spawn
                                        cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start( cluster.domains[d].cuda_pinned_buff.data(), cluster.domains[d].cuda_pinned_buff.data(),'N',0 );
				} else {
					//cilk_spawn
                                        cluster.domains[d].B1Kplus.DenseMatVecCUDA_wo_Copy_start_fl( cluster.domains[d].cuda_pinned_buff_fl.data(), cluster.domains[d].cuda_pinned_buff_fl.data(),'N',0 );
				}
		  	}
		}

		#pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++) {
		  	if (!configuration.combine_sc_and_spds) {
				if (cluster.domains[d].B1Kplus.is_on_acc == 0) {
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


void IterSolverGPU::Apply_Prec( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ) {

	 time_eval.totalTime.start();

	 time_eval.timeEvents[0].start();

	#pragma omp parallel for
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		
		// TODO_GPU: Zde se tvori vektor za danou domenu, kterym se nasobi matice predpodminovace - tohle nechame na CPU. Vektor x_in je za cely cluster (vice domen, ktere sdili pamet na uzlu) 
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d]->B1_comp_dom.rows, 0.0 );
		for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_in[ cluster.domains[d]->lambda_map_sub_local[i]] * cluster.domains[d]->B1_scale_vec[i]; // includes B1 scaling

		switch (USE_PREC) {
		case FETIConfiguration::PRECONDITIONER::LUMPED:
			cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'T');
			cluster.domains[d]->K.MatVec(*cluster.x_prim_cluster1[d], *cluster.x_prim_cluster2[d],'N');
			if (cluster.domains[d]->_RegMat.nnz > 0) {
				cluster.domains[d]->_RegMat.MatVecCOO(*cluster.x_prim_cluster1[d], *cluster.x_prim_cluster2[d],'N', 1.0, -1.0);
			}
			break;
		case FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION:
			cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster2[d], 'T');
			break;
		case FETIConfiguration::PRECONDITIONER::DIRICHLET:
			// TODO_GPU
			// zde je verze predpominovace na CPU, ktera se ma predelat pro GPU. Pro zacatek udelame na GPU pouze DIRICHLET
			// nasledujici radek by mel zustat na CPU 
			if(cluster.domains[d]->Prec.is_on_acc == 1) {
					cluster.domains[d]->B1t_DirPr.MatVec (x_in_tmp, cluster.domains[d]->cuda_pinned_buff, 'N');
			} else {
					cluster.domains[d]->B1t_DirPr.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'N');
			}
			// nasledujici readek by mel jit na GPU - je otazkou jak udelal OpenMP paralelizaci vs. ovladani CUDY a streamu - podivat se na reseni v Apply_A a pro zacatek obe reseni sjednotit 
//			cluster.domains[d]->Prec.DenseMatVec(*cluster.x_prim_cluster1[d], *cluster.x_prim_cluster2[d],'N');
			break;
		case FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET:
			cluster.domains[d]->B1t_DirPr.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'N');
			cluster.domains[d]->Prec.MatVec(*cluster.x_prim_cluster1[d], *cluster.x_prim_cluster2[d],'N');
			break;
		case FETIConfiguration::PRECONDITIONER::MAGIC:
			cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'T');
			cluster.domains[d]->Prec.MatVec(*cluster.x_prim_cluster1[d], *cluster.x_prim_cluster2[d],'N');
			break;
		case FETIConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::error("Not implemented preconditioner.\n");
		}

	}

	// Only for DIRICHLET case for now
	if(USE_PREC == FETIConfiguration::PRECONDITIONER::DIRICHLET) {
		PUSH_RANGE("Apply_Prec GEMV", 3)
		// For those domains that did fit in to GPU initialize Prec application
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if(cluster.domains[d]->Prec.is_on_acc == 1) {
                cluster.domains[d]->Prec.DenseMatVecCUDA_wo_Copy_start( cluster.domains[d]->cuda_pinned_buff.data(), cluster.domains[d]->cuda_pinned_buff.data(),'N',0 );
			}
		}
		// For those domains that didn't fit in to GPU apply Prec on CPU
		#pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++)
		{
			if(cluster.domains[d]->Prec.is_on_acc == 0) {
				cluster.domains[d]->Prec.DenseMatVec(*cluster.x_prim_cluster1[d], *cluster.x_prim_cluster2[d],'N');
			}
		}
		// Wait for GPU to finish operations
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (cluster.domains[d]->Prec.is_on_acc == 1) {
                cluster.domains[d]->Prec.DenseMatVecCUDA_wo_Copy_sync ( );
			}
		}
		POP_RANGE // Apply_Prec
	}

	std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
	SEQ_VECTOR < double > y_out_tmp;
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		y_out_tmp.resize( cluster.domains[d]->B1_comp_dom.rows );

		switch (USE_PREC) {
		case FETIConfiguration::PRECONDITIONER::LUMPED:
		case FETIConfiguration::PRECONDITIONER::WEIGHT_FUNCTION:
		case FETIConfiguration::PRECONDITIONER::MAGIC:
			cluster.domains[d]->B1_comp_dom.MatVec (*cluster.x_prim_cluster2[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
		//TODO  check if MatVec is correct (DenseMatVec!!!)
		case FETIConfiguration::PRECONDITIONER::DIRICHLET:
			if(cluster.domains[d]->Prec.is_on_acc == 1) {
				cluster.domains[d]->B1t_DirPr.MatVec (cluster.domains[d]->cuda_pinned_buff, y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			} else {
				cluster.domains[d]->B1t_DirPr.MatVec (*cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			}
			break;
		case FETIConfiguration::PRECONDITIONER::SUPER_DIRICHLET:
			cluster.domains[d]->B1t_DirPr.MatVec (*cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
		case FETIConfiguration::PRECONDITIONER::NONE:
			break;
		default:
			eslog::error("Not implemented preconditioner.\n");
		}

		for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++)
			cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += y_out_tmp[i] * cluster.domains[d]->B1_scale_vec[i]; // includes B1 scaling
	}
	 time_eval.timeEvents[0].end();

	 time_eval.timeEvents[1].start();
	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
	 time_eval.timeEvents[1].end();

	 time_eval.totalTime.end();
}


void IterSolverGPU::apply_A_l_comp_dom_B_P( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {

    	eslog::error("apply_A_l_comp_dom_B_P - is not implemented in itersolverGPU.cpp .\n");

		// See cpu/itersolver.cpu for reference implementation 

}

void IterSolverGPU::apply_A_l_comp_dom_B_P_local( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {

    	eslog::error("apply_A_l_comp_dom_B_P_local - is not implemented in itersolverGPU.cpp .\n");

		// See cpu/itersolver.cpu for reference implementation 

}

void IterSolverGPU::apply_A_l_comp_dom_B_P_local_sparse( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<esint> & in_indices, SEQ_VECTOR<double> & in_values, SEQ_VECTOR<esint> & out_indices, SEQ_VECTOR<double> & out_values) {
//	 time_eval.totalTime.start();
    if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {

    	eslog::error("Hybrid FETI not supported in apply_A_l_comp_dom_B_P_local_sparse in itersolverGPU.cpp.\n");

		// See cpu/itersolver.cpu for reference implementation 

    }

    if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 0) {
//    	 time_eval.timeEvents[0].start();
    	SEQ_VECTOR <int> is_empty ( cluster.x_prim_cluster1.size(), true);
                #pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR <double> x_in_tmp( cluster.domains[d]->B1_comp_dom.rows, 0.0 );

			esint iter_in = 0;
			esint iter_lm = 0;

			do {

				if (cluster.domains[d]->lambda_map_sub[iter_lm] == in_indices[iter_in]) {
					x_in_tmp[iter_lm] = in_values[iter_in];
					iter_in++;
					iter_lm++;
					is_empty[d] = false;
				} else {
					if (cluster.domains[d]->lambda_map_sub[iter_lm] > in_indices[iter_in]) {
						iter_in++;
					} else {
						iter_lm++;
					}
				}

			} while (iter_lm != (esint)cluster.domains[d]->lambda_map_sub.size() && iter_in != (esint)in_indices.size() );

			if (!is_empty[d]) {
				//CPU
				// cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'T');
				//GPU - LSC 
				if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
					for (size_t i = 0; i < x_in_tmp.size(); i++) {
						cluster.domains[d]->cuda_pinned_buff[i] = x_in_tmp[i];
					}
				} else {
					for (size_t i = 0; i < x_in_tmp.size(); i++) {
						cluster.domains[d]->compressed_tmp[i] = x_in_tmp[i];
					}
				}
			}
		}


//         time_eval.timeEvents[0].end();

//         time_eval.timeEvents[1].start();
        if (cluster.USE_HFETI == 0) {
        	//cluster.multKplusFETI (cluster.x_prim_cluster1);
			PUSH_RANGE("apply_A_l_comp_dom_B_P_local_sparse GEMV", 4)
			#pragma omp parallel for
			for (size_t d = 0; d < cluster.domains.size(); d++) {
				if (!is_empty[d] && (cluster.domains[d]->B1Kplus.is_on_acc == 1)) {

					// Each thread needs to call this otherwise non-master threads will choose default GPU
					// performance consideration: swap with CuCtxSetCurrent() call ~ 2x faster
					cudaSetDevice(cluster.domains[d]->B1Kplus.device_id);
				// CPU
				// cluster.domains[d]->multKplusLocal(*cluster.x_prim_cluster1[d]);
				// GPU - LSC
#ifdef SHARE_SC
					cluster.domains[d]->B1Kplus.DenseMatVecCUDA_shared_wo_Copy_start(cluster.domains[d]->cuda_pinned_buff, cluster.domains[d]->cuda_pinned_buff,'N',0);
#else
					cluster.domains[d]->B1Kplus.DenseMatVecCUDA_wo_Copy_start(cluster.domains[d]->cuda_pinned_buff.data(), cluster.domains[d]->cuda_pinned_buff.data(),'N',0);
#endif
				}
			}
		} else {
			eslog::error("apply_A_l_comp_dom_B_P_local_sparse - ERROR HTFETI_1x .\n");
			// cluster.multKplusHFETI(cluster.x_prim_cluster1);
        }

	#pragma omp parallel for
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		// If non-empty and didnt fit into GPU, compute on CPU. TODO check if we can use *cluster.x_prim_cluster1[d]!!!
		if ((!is_empty[d]) && (cluster.domains[d]->B1Kplus.is_on_acc == 0) ) {
			if (configuration.combine_sc_and_spds)
			{
				cluster.domains[d]->B1_comp_dom.MatVec (cluster.domains[d]->compressed_tmp, *cluster.x_prim_cluster1[d], 'T');
				cluster.domains[d]->multKplusLocal(*cluster.x_prim_cluster1[d]);
				cluster.domains[d]->B1_comp_dom.MatVec (*cluster.x_prim_cluster1[d], cluster.domains[d]->compressed_tmp2, 'N', 0, 0, 0.0);
			}
			else
			{
				cluster.domains[d]->B1Kplus.DenseMatVec ( cluster.domains[d]->compressed_tmp, cluster.domains[d]->compressed_tmp2);
			}
		}
	}
//         time_eval.timeEvents[1].end();

//         time_eval.timeEvents[2].start();
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);

		SEQ_VECTOR < double > y_out_tmp;
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (!is_empty[d]) {
				// CPU
				// y_out_tmp.resize( cluster.domains[d]->B1_comp_dom.rows );
				// cluster.domains[d]->B1_comp_dom.MatVec (*cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
				// for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++)
					// cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += y_out_tmp[i];

				// GPU - LSC 
				if (cluster.domains[d]->B1Kplus.is_on_acc == 1) {
					cluster.domains[d]->B1Kplus.DenseMatVecCUDA_wo_Copy_sync ( );
					for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += cluster.domains[d]->cuda_pinned_buff[i];					
					}
				} else {
					for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++) {
						cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += cluster.domains[d]->compressed_tmp2[i];
					}
				}
			}
		}
		POP_RANGE // apply_A_l_comp_dom_B_P_local_sparse

//		time_eval.timeEvents[2].end();
    }


    if (cluster.USE_KINV == 0) {
//    	 time_eval.timeEvents[0].start();
    	SEQ_VECTOR <int> is_empty ( cluster.x_prim_cluster1.size(), true);
 		#pragma omp parallel for
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d]->B1_comp_dom.rows, 0.0 );

			esint iter_in = 0;
			esint iter_lm = 0;

			do {

				if (cluster.domains[d]->lambda_map_sub[iter_lm] == in_indices[iter_in]) {
					x_in_tmp[iter_lm] = in_values[iter_in];
					iter_in++;
					iter_lm++;
					is_empty[d] = false;
				} else {
					if (cluster.domains[d]->lambda_map_sub[iter_lm] > in_indices[iter_in]) {
						iter_in++;
					} else {
						iter_lm++;
					}
				}

			} while (iter_lm != (esint)cluster.domains[d]->lambda_map_sub.size() && iter_in != (esint)in_indices.size() );

			if (!is_empty[d])
				cluster.domains[d]->B1_comp_dom.MatVec (x_in_tmp, *cluster.x_prim_cluster1[d], 'T');
		}

//         time_eval.timeEvents[0].end();

//         time_eval.timeEvents[1].start();
        if (cluster.USE_HFETI == 0) {
        	//cluster.multKplusFETI (cluster.x_prim_cluster1);

			#pragma omp parallel for
			for (size_t d = 0; d < cluster.domains.size(); d++) {
				if (!is_empty[d]) {
					cluster.domains[d]->multKplusLocal(*cluster.x_prim_cluster1[d]);
				}
			}
		} else {
    		cluster.multKplusHFETI(cluster.x_prim_cluster1);
        }
//         time_eval.timeEvents[1].end();


//         time_eval.timeEvents[2].start();
		std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);

		SEQ_VECTOR < double > y_out_tmp;
		for (size_t d = 0; d < cluster.domains.size(); d++) {
			if (!is_empty[d]) {
				y_out_tmp.resize( cluster.domains[d]->B1_comp_dom.rows );
				cluster.domains[d]->B1_comp_dom.MatVec (*cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

				for (size_t i = 0; i < cluster.domains[d]->lambda_map_sub_local.size(); i++)
					cluster.compressed_tmp[ cluster.domains[d]->lambda_map_sub_local[i] ] += y_out_tmp[i];
			}
		}

//		time_eval.timeEvents[2].end();
    }

	out_indices = cluster.my_lamdas_indices;
	out_values  = cluster.compressed_tmp;

//     time_eval.totalTime.end();

}
