/*
 * clusterGPU.h
 *
 *  Created on: Feb 24, 2016
 *      Author: lriha
 */

#ifndef SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_


#include "../cluster.h"
#include <cuda.h>
#include <cuda_runtime.h>
namespace espreso {

class ClusterGPU: public ClusterBase
{

public:
	// Constructor
	ClusterGPU(const ESPRESOSolver &configuration, Instance *instance_in): ClusterBase(configuration, instance_in) { }
	~ClusterGPU();

	void Create_SC_perDomain( bool USE_FLOAT );
	void GetSchurComplement( bool USE_FLOAT, eslocal i );

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

#ifdef SHARE_SC
	SEQ_VECTOR <double *> SC_dense_val_orig;
	SEQ_VECTOR <float *> SC_dense_val_orig_fl;
#endif
};

}



#endif /* SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_ */
