/*
 * clusterGPU.h
 *
 *  Created on: Feb 24, 2016
 *      Author: lriha
 */

#ifndef SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_


#include "feti/specific/cluster.h"
#include <cuda.h>
#include <cuda_runtime.h>
namespace espreso {

class ClusterGPU: public ClusterBase
{

public:
	// Constructor
	ClusterGPU(const FETIConfiguration &configuration, DataHolder *instance_in): ClusterBase(configuration, instance_in), device_id(-1) { }
	~ClusterGPU();

	void Create_SC_perDomain( bool USE_FLOAT );
        void CreateDirichletPrec( DataHolder *instance );

	void GetSchurComplement( bool USE_FLOAT, esint i );
        void GetDirichletPrec( DataHolder *instance, esint d );


	void SetupKsolvers ( );

	void multKplusGlobal_GPU   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );


	void multKplus_HF      (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out);
	void multKplus_HF_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);

	void multKplus_HF_Loop1 (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_CP    ();

	void multKplus_HF_Loop2_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out);
	void multKplus_HF_Loop2_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_Loop2_MIX  (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);

	void CreateCudaStreamPool();
	void DestroyCudaStreamPool();

	SEQ_VECTOR <cudaStream_t> cuda_stream_pool;
private:
	void GetGPU();

        size_t  GPU_free_mem;
        size_t  GPU_total_mem;

	int device_id;



#ifdef SHARE_SC
	SEQ_VECTOR <double *> SC_dense_val_orig;
	SEQ_VECTOR <float *> SC_dense_val_orig_fl;
#endif
};

}



#endif /* SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_ */
