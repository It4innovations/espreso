
#include "itersolveracc.h"


using namespace espreso;
// *** Action of K+ routines *********************************************

void IterSolverAcc::apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out) {
    time_eval.totalTime.start();
    // number of Xeon Phi devices 
    eslocal numDevices = config::solver::N_MICS;

    if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 1) {
        // HFETI on MIC using Schur

        time_eval.timeEvents[0].start();

        eslocal maxDevNumber = numDevices;

        // *** Part 1.1 - prepare vectors for FETI operator with SC
        // at first - work with domains assigned to MICs
        for ( eslocal i = 0; i < maxDevNumber; i++ ) {
            cilk_for (eslocal d = 0; d < cluster.accDomains[i].size(); ++d) {
                eslocal domN = cluster.accDomains[i].at(d);
                for ( eslocal j = 0; j < cluster.domains[domN].lambda_map_sub_local.size(); j++ ) {
                    cluster.B1KplusPacks[i].SetX(d, j, x_in[ cluster.domains[domN].lambda_map_sub_local[j]]);
                    cluster.domains[domN].compressed_tmp2[j] = x_in[ cluster.domains[domN].lambda_map_sub_local[j]];
                }
                // *** Part 2.1 - prepare vectors for HFETI operator on CPU
                cluster.domains[domN].B1_comp_dom.MatVec (cluster.domains[domN].compressed_tmp2, cluster.x_prim_cluster1[domN], 'T');
            }
        }

        // *** Part 1.2 work with domains staying on CPU
        cilk_for (eslocal d = 0; d < cluster.hostDomains.size(); ++d ) {
            eslocal domN = cluster.hostDomains.at(d);
            for ( eslocal j = 0; j < cluster.domains[domN].lambda_map_sub_local.size(); j++ ) {
                cluster.domains[domN].compressed_tmp2[j] = x_in[ cluster.domains[domN].lambda_map_sub_local[j]];
            }
            // *** Part 2.2 - prepare vectors for HFETI operator on CPU
            cluster.domains[domN].B1_comp_dom.MatVec (cluster.domains[domN].compressed_tmp2, cluster.x_prim_cluster1[domN], 'T');
        }

        // *** Part 3 - execute FETI SC operator 
        // spawn the computation on MICs
#pragma omp parallel num_threads( maxDevNumber )
        {
            // start async. computation on MICs
            if (cluster.accDomains[omp_get_thread_num()].size() > 0) {
                cluster.B1KplusPacks[ omp_get_thread_num() ].DenseMatsVecsMIC_Start( 'N' );
            }
        }
        // meanwhile compute the same for domains staying on CPU
        double startCPU = Measure::time();
        cilk_for (eslocal d = 0; d < cluster.hostDomains.size(); ++d ) {
            eslocal domN = cluster.hostDomains.at(d);
            cluster.domains[domN].B1Kplus.DenseMatVec( cluster.domains[domN].compressed_tmp2, cluster.domains[domN].compressed_tmp);
        }

        for ( eslocal i = 0 ; i < numDevices; ++i ) {
            cluster.B1KplusPacks[ i ].DenseMatsVecsRestCPU( 'N' );    
            long start = (long) (cluster.B1KplusPacks[i].getNMatrices()*cluster.B1KplusPacks[i].getMICratio());
            cilk_for (  long d = start ; d < cluster.B1KplusPacks[i].getNMatrices(); ++d ) {
                cluster.B1KplusPacks[i].GetY(d, cluster.domains[cluster.accDomains[i].at(d)].compressed_tmp);
            }

        }

        time_eval.timeEvents[0].end();

        // *** Part 4 - Execute HTFETI operator
        // perform simultaneous computation on CPU
        time_eval.timeEvents[1].start();
        cluster.multKplusGlobal_Kinv( cluster.x_prim_cluster1 );
        time_eval.timeEvents[1].end();
        time_eval.timeEvents[2].start();
        double CPUtime = Measure::time() - startCPU;


        // *** Part 5 - Finalize transfers of the result of the FETI SC operator
        // from MICs back to CPU
#pragma omp parallel num_threads( maxDevNumber )
        {
            // synchronize computation
            if (cluster.accDomains[omp_get_thread_num()].size() > 0) {
                cluster.B1KplusPacks[ omp_get_thread_num() ].DenseMatsVecsMIC_Sync(  );
            }
        }

        // extract the result from MICs
        for ( eslocal i = 0; i < maxDevNumber; i++ ) {
            long end = (long) (cluster.B1KplusPacks[i].getNMatrices()*cluster.B1KplusPacks[i].getMICratio()); 
            cilk_for ( eslocal d = 0 ; d < end; ++d ) {
                cluster.B1KplusPacks[i].GetY(d, cluster.domains[cluster.accDomains[i].at(d)].compressed_tmp);
            }
        }
        // *** Part 6 - Lambda values, results of the FETI SC operator, are
        // combined back into single lambda vector per cluster -
        // cluster.compressed_tmp
        std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
        for (eslocal d = 0; d < cluster.domains.size(); ++d) {
            for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
        }

        // *** Part 7 - calculate lambda values from the results of the HTFETI
        // operator (results are in cluster.x_prim_cluster1)
        //      //              using cluster.domains[d].B1_comp_dom gluing
        //      matrix and store them in the temp. vector y_out_tmp.
        //              //  
        SEQ_VECTOR < double > y_out_tmp;
        for (eslocal d = 0; d < cluster.domains.size(); d++) {
            y_out_tmp.resize( cluster.domains[d].B1_comp_dom.rows );
            cluster.domains[d].B1_comp_dom.MatVec (cluster.x_prim_cluster1[d], y_out_tmp, 'N', 0, 0, 0.0); // will add (summation per elements) all partial results into y_out

            for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += y_out_tmp[i];
        }

        if ( config::solver::LOAD_BALANCING ) {
            // update the ratio between the cpu and mic
            double r = cluster.B1KplusPacks[0].getMICratio();
            double MICtime = cluster.B1KplusPacks[0].getElapsedTime();
            double newRatio = (r * CPUtime) / (r * CPUtime + MICtime * (1 - r));
            //std::cout << "TEST " << r << " " <<  CPUtime<< " "  << MICtime << " " << newRatio << std::endl;

#pragma omp parallel num_threads( maxDevNumber )
            {
                cluster.B1KplusPacks[omp_get_thread_num()].setMICratio( newRatio );
            }
        }

        time_eval.timeEvents[2].end();
    }

    if (cluster.USE_KINV == 1 && cluster.USE_HFETI == 0) {
        // classical FETI on MIC using Schur

        time_eval.timeEvents[0].start();
        time_eval.timeEvents[0].end();

        time_eval.timeEvents[1].start();

        eslocal maxDevNumber = numDevices;

        // *** Part 1.1 - prepare vectors for FETI operator with SC
        // at first - work with domains assigned to MICs
        for ( eslocal i = 0; i < maxDevNumber; i++ ) {
            cilk_for (eslocal d = 0; d < cluster.accDomains[i].size(); d++) {
                eslocal domN = cluster.accDomains[i].at(d);
                for ( eslocal j = 0; j < cluster.domains[domN].lambda_map_sub_local.size(); j++ ) {
                    cluster.B1KplusPacks[i].SetX(d, j, x_in[ cluster.domains[domN].lambda_map_sub_local[j]]);
            }
            }
        }
        
        // *** Part 1.2 work with domains staying on CPU
        cilk_for (eslocal d = 0; d < cluster.hostDomains.size(); ++d ) {
            eslocal domN = cluster.hostDomains.at(d);
            for ( eslocal j = 0; j < cluster.domains[domN].lambda_map_sub_local.size(); j++ ) {
                cluster.domains[domN].compressed_tmp2[j] = x_in[ cluster.domains[domN].lambda_map_sub_local[j]];
            }
        }

//#pragma omp parallel num_threads( maxDevNumber )
//        {
//            // start async. computation on MICs
for (int mic=0; mic<maxDevNumber; mic++)            
            if (cluster.accDomains[mic].size() > 0) {
                cluster.B1KplusPacks[ mic ].DenseMatsVecsMIC_Start( 'N' );
            }
//        }
        // meanwhile compute the same for domains staying on CPU
        double startCPU = Measure::time();
        cilk_for (eslocal d = 0; d < cluster.hostDomains.size(); ++d ) {
            eslocal domN = cluster.hostDomains.at(d);
            cluster.domains[domN].B1Kplus.DenseMatVec( cluster.domains[domN].compressed_tmp2, cluster.domains[domN].compressed_tmp);
        }
        for ( eslocal i = 0 ; i < numDevices; ++i ) {
            cluster.B1KplusPacks[ i ].DenseMatsVecsRestCPU( 'N' );    
            long start = (long) (cluster.B1KplusPacks[i].getNMatrices()*cluster.B1KplusPacks[i].getMICratio());
            cilk_for (  long d = start; d < cluster.B1KplusPacks[i].getNMatrices(); ++d ) {
                cluster.B1KplusPacks[i].GetY(d, cluster.domains[cluster.accDomains[i].at(d)].compressed_tmp);
            }
        }

        double CPUtime = Measure::time() - startCPU;

        // *** Part 5 - Finalize transfers of the result of the FETI SC operator
        // from MICs back to CPU
//#pragma omp parallel num_threads( maxDevNumber )
//        {
            // synchronize computation
            for (int mic=0 ; mic < maxDevNumber; mic++)
            if (cluster.accDomains[mic].size() > 0) {
                cluster.B1KplusPacks[ mic ].DenseMatsVecsMIC_Sync(  );
            }
//        }

        // extract the result from MICs
        for ( eslocal i = 0; i < maxDevNumber; i++ ) {
            long end = (long) (cluster.B1KplusPacks[i].getNMatrices()*cluster.B1KplusPacks[i].getMICratio()); 
            cilk_for ( eslocal d = 0 ; d < end; ++d ) {
                cluster.B1KplusPacks[i].GetY(d, cluster.domains[cluster.accDomains[i].at(d)].compressed_tmp);
            }
        }
        time_eval.timeEvents[1].end();

        // *** Part 6 - Lambda values, results of the FETI SC operator, are
        // combined back into single lambda vector per cluster -
        // cluster.compressed_tmp
        time_eval.timeEvents[2].start();
        std::fill( cluster.compressed_tmp.begin(), cluster.compressed_tmp.end(), 0.0);
        for (eslocal d = 0; d < cluster.domains.size(); ++d) {
            for (eslocal i = 0; i < cluster.domains[d].lambda_map_sub_local.size(); i++)
                cluster.compressed_tmp[ cluster.domains[d].lambda_map_sub_local[i] ] += cluster.domains[d].compressed_tmp[i];
        }
       
       if ( config::solver::LOAD_BALANCING ) {
            // update the ratio between the cpu and mic
            double r = cluster.B1KplusPacks[0].getMICratio();
            double MICtime = cluster.B1KplusPacks[0].getElapsedTime();
            double newRatio = (r * CPUtime) / (r * CPUtime + MICtime * (1 - r));
            std::cout << "TEST " << r << " " <<  CPUtime<< " "  << MICtime << " " << newRatio << std::endl;

#pragma omp parallel num_threads( maxDevNumber )
            {
                cluster.B1KplusPacks[omp_get_thread_num()].setMICratio( newRatio );
            }
        }


        time_eval.timeEvents[2].end();
        
    }

/*    if (cluster.USE_KINV == 0) {
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
            vectorsPerAcc.reserve(config::solver::N_MICS);
            SEQ_VECTOR<eslocal> nVecPerMIC;
            nVecPerMIC.resize(config::solver::N_MICS);

            for (eslocal i = 0; i < config::solver::N_MICS; i++) {
                nVecPerMIC[i] = nMatrices / config::solver::N_MICS;
            }

            for (eslocal i = 0 ; i < nMatrices % config::solver::N_MICS; i++ ) {
                nVecPerMIC[i]++;
            }

            eslocal offset = 0;
            for (eslocal i = 0; i < config::solver::N_MICS; i++) {
                vectorsPerAcc[i] = new SEQ_VECTOR<double>*[ nVecPerMIC[ i ] ];
                for (eslocal j = offset; j < offset + nVecPerMIC[ i ]; j++) {
                    (vectorsPerAcc[i])[j - offset] = &(cluster.x_prim_cluster1[j]); 
                }
                offset += nVecPerMIC[i];
                //     cluster.solver[i].Solve(vectorsPerAcc[i]);
            }

#pragma omp parallel num_threads(config::solver::N_MICS)
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

    }*/

    time_eval.timeEvents[3].start();
    All_Reduce_lambdas_compB(cluster, cluster.compressed_tmp, y_out);
    time_eval.timeEvents[3].end();

    time_eval.totalTime.end();

}


