#ifndef SOLVER_SPECIFIC_ACC_ITERSOLVERGPU_H_
#define SOLVER_SPECIFIC_ACC_ITERSOLVERGPU_H_

#include "../itersolver.h"

namespace espreso {

class IterSolverGPU: public IterSolverBase
{
	   virtual void apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double>  & y_out);

	   void apply_A_FETI_SC_forHFETI   ( Cluster & cluster );

	   void apply_A_FETI_SC_forFETI   ( Cluster & cluster, SEQ_VECTOR<double> & x_in );

	   void apply_A_htfeti_SC ( Cluster & cluster );



};

}


#endif /* SOLVER_SPECIFIC_ACC_ITERSOLVERGPU_H_ */
