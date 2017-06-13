
#include "itersolvercpu.h"

using namespace espreso;

// *** Action of K+ routines *********************************************

void IterSolverCPU::apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {
       time_eval.totalTime.start();

    if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {
         time_eval.timeEvents[0].start();
        #pragma omp parallel for
        for (size_t d = 0; d < cluster.domains.size(); d++) {
            //SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
            for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++) {
                cluster.domains[d].compressed_tmp2[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
            }

            cluster.domains[d].B1Kplus.DenseMatVec (cluster.domains[d].compressed_tmp2, cluster.domains[d].compressed_tmp);

            cluster.domains[d].B1_comp_dom.MatVec (cluster.domains[d].compressed_tmp2, cluster.x_prim_cluster1[d], 'T');
        }
         time_eval.timeEvents[0].end();

         time_eval.timeEvents[1].start();
        cluster.multKplusGlobal_Kinv( cluster.x_prim_cluster1 );
         time_eval.timeEvents[1].end();

         time_eval.timeEvents[2].start();
        std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);

        for (size_t d = 0; d < cluster.domains.size(); d++) {
            for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
        }

        //std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
        SEQ_VECTOR < double > y_out_tmp;
        for (size_t d = 0; d < cluster.domains.size(); d++) {
            y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
            cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

            for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
        }
         time_eval.timeEvents[2].end();

    }


    if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 0) {
         time_eval.timeEvents[0].start();
        //// POZOR - jen pro porovnani vykonu CPU a GPU
        //cilk_for (size_t d = 0; d < cluster.domains.size(); d++) {
        //	SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
        //	for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
        //		x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
        //	cluster.domains[d].B1Kplus.DenseMatVec(x_in_tmp, cluster.domains[d].compressed_tmp);
        //}
        //// POZOR - jen pro porovnani vykonu CPU a GPU
         time_eval.timeEvents[0].end();

         time_eval.timeEvents[1].start();
        #pragma omp parallel for
        for (size_t d = 0; d < cluster.domains.size(); d++) {
            SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows );
            for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
            cluster.domains[d].B1Kplus.DenseMatVec            (x_in_tmp, cluster.domains[d].compressed_tmp);
        }
         time_eval.timeEvents[1].end();

         time_eval.timeEvents[2].start();
        std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
        for (size_t d = 0; d < cluster.domains.size(); d++) {
            for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
        }
         time_eval.timeEvents[2].end();
    }


    if (cluster.USE_KINV == 0) {
         time_eval.timeEvents[0].start();
        #pragma omp parallel for
        for (size_t d = 0; d < cluster.domains.size(); d++) {
            SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
            for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]];
            cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
            //cluster.x_prim_cluster2[d] = cluster.x_prim_cluster1[d]; // POZOR zbytecne kopirovani // prim norm
        }
         time_eval.timeEvents[0].end();

         time_eval.timeEvents[1].start();
        if (cluster.USE_HFETI == 0) {
            #pragma omp parallel for
        	for (size_t d = 0; d < cluster.domains.size(); d++) {
                cluster.domains[d].multKplusLocal(cluster.x_prim_cluster1[d]);
            }
		} else {
			cluster.multKplusGlobal_l(cluster.x_prim_cluster1);
        }
         time_eval.timeEvents[1].end();


         time_eval.timeEvents[2].start();
        std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
        SEQ_VECTOR < double > y_out_tmp;
        for (size_t d = 0; d < cluster.domains.size(); d++) {
            y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
            cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

            for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
        }
         time_eval.timeEvents[2].end();

    }

     time_eval.timeEvents[3].start();
    All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
     time_eval.timeEvents[3].end();

     time_eval.totalTime.end();

}

void IterSolverCPU::apply_prec_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out ) {

	time_eval.totalTime.start();

	time_eval.timeEvents[0].start();

	#pragma omp parallel for
	for (size_t d = 0; d < cluster.domains.size(); d++) {
		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
			x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling

		switch (USE_PREC) {
		case ESPRESO_PRECONDITIONER::LUMPED:
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
			cluster.domains[d].K.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			if (cluster.domains[d]._RegMat.nnz > 0) {
				cluster.domains[d]._RegMat.MatVecCOO(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N', 1.0, -1.0);
			}
			break;
		case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster2[d], 'T');
			break;
		case ESPRESO_PRECONDITIONER::DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
			//cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			cluster.domains[d].Prec.DenseMatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			break;
		case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
//			cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			break;
		case ESPRESO_PRECONDITIONER::MAGIC:
			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
			cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
			break;
		case ESPRESO_PRECONDITIONER::NONE:
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
		case ESPRESO_PRECONDITIONER::LUMPED:
		case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
		case ESPRESO_PRECONDITIONER::MAGIC:
			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
    //TODO  check if MatVec is correct (DenseMatVec!!!) 
		case ESPRESO_PRECONDITIONER::DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
		case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
			cluster.domains[d].B1t_DirPr.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
			break;
		case ESPRESO_PRECONDITIONER::NONE:
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


//	time_eval.totalTime.start();
//
//	time_eval.timeEvents[0].start();
//
//	#pragma omp parallel for
//	for (size_t d = 0; d < cluster.domains.size(); d++) {
//		SEQ_VECTOR < double > x_in_tmp ( cluster.domains[d].B1_comp_dom.rows, 0.0 );
//		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
//			x_in_tmp[i] = x_in[ cluster.domains[d].lambda_map_sub_local[i]] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling
//
//		switch (USE_PREC) {
//		case ESPRESO_PRECONDITIONER::LUMPED:
//			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
//			cluster.domains[d].K.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
//			cluster.domains[d]._RegMat.MatVecCOO(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N', -1.0);
//			break;
//		case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
//			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster2[d], 'T');
//			break;
//		case ESPRESO_PRECONDITIONER::DIRICHLET:
//			cluster.domains[d].B1t_DirPr.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
//			//cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
//			cluster.domains[d].Prec.DenseMatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
//			break;
//		case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
//			cluster.domains[d].B1t_DirPr.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'N');
//			cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
//			break;
//		case ESPRESO_PRECONDITIONER::MAGIC:
//			cluster.domains[d].B1_comp_dom.MatVec (x_in_tmp, cluster.x_prim_cluster1[d], 'T');
//			cluster.domains[d].Prec.MatVec(cluster.x_prim_cluster1[d], cluster.x_prim_cluster2[d],'N');
//			break;
//		case ESPRESO_PRECONDITIONER::NONE:
//			break;
//		default:
//			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
//		}
//
//	}
//
//	std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
//	SEQ_VECTOR < double > y_out_tmp;
//	for (size_t d = 0; d < cluster.domains.size(); d++) {
//		y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
//
//
//		switch (USE_PREC) {
//		case ESPRESO_PRECONDITIONER::LUMPED:
//		case ESPRESO_PRECONDITIONER::WEIGHT_FUNCTION:
//		case ESPRESO_PRECONDITIONER::MAGIC:
//			cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
//			break;
//    //TODO  check if MatVec is correct (DenseMatVec!!!)
//		case ESPRESO_PRECONDITIONER::DIRICHLET:
//			cluster.domains[d].B1t_DirPr.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
//			break;
//		case ESPRESO_PRECONDITIONER::SUPER_DIRICHLET:
//			cluster.domains[d].B1t_DirPr.MatVec (cluster.x_prim_cluster2[d], y_out_tmp, 'T', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out
//			break;
//		case ESPRESO_PRECONDITIONER::NONE:
//			break;
//		default:
//			ESINFO(GLOBAL_ERROR) << "Not implemented preconditioner.";
//		}
//
//
//		for (size_t i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
//			cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i] * cluster.domains[d].B1_scale_vec[i]; // includes B1 scaling
//	}
//	time_eval.timeEvents[0].end();
//
//
//	time_eval.timeEvents[1].start();
//	All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
//	time_eval.timeEvents[1].end();
//
//
//	time_eval.totalTime.end();

}



