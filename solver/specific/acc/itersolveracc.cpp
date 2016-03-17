
#include "itersolveracc.h"



// *** Action of K+ routines *********************************************

void IterSolverAcc::apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {
    time_eval.totalTime.start();

    // number of Xeon Phi devices (0 for fallback to CPU)
    eslocal numDevices = N_MICS;


    if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {

std::cout << "INITIATING HFETI"<<std::endl;
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
            //cilk_for (eslocal d = 0; d < cluster.domains.size(); d++) {
            eslocal nMatrices = cluster.domains.size(); 
            SEQ_VECTOR<SEQ_VECTOR<double>**> vectorsPerAcc;
            vectorsPerAcc.reserve(N_MICS);
            SEQ_VECTOR<eslocal> nVecPerMIC;
            nVecPerMIC.resize(N_MICS);

            for (eslocal i = 0; i < N_MICS; i++) {
                nVecPerMIC[i] = nMatrices / N_MICS;
            }

            for (eslocal i = 0 ; i < nMatrices % N_MICS; i++ ) {
                nVecPerMIC[i]++;
            }

            eslocal offset = 0;
            for (eslocal i = 0; i < N_MICS; i++) {
                vectorsPerAcc[i] = new SEQ_VECTOR<double>*[ nVecPerMIC[ i ] ];
                for (eslocal j = offset; j < offset + nVecPerMIC[ i ]; j++) {
                    (vectorsPerAcc[i])[j - offset] = &(cluster.x_prim_cluster1[j]); 
                }
                offset += nVecPerMIC[i];
           //     cluster.solver[i].Solve(vectorsPerAcc[i]);
            }

#pragma omp parallel num_threads(N_MICS)
            {
                eslocal myAcc = omp_get_thread_num();
                cluster.solver[myAcc].Solve(vectorsPerAcc[myAcc]);
                delete [] vectorsPerAcc[myAcc];
            }
            //for(eslocal i = 0; i < cluster.domains.size(); i++)
            //cluster.domains[i].multKplusLocal(cluster.x_prim_cluster1[i]);
            // }
        } else {
            // TODO NOT YET IMPLEMENTED ON MIC
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


