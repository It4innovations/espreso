
#ifndef SOLVER_SPECIFIC_ACC_ITERSOLVERACC_H_
#define SOLVER_SPECIFIC_ACC_ITERSOLVERACC_H_


#include "../itersolver.h"

class IterSolverAcc: public IterSolverBase
{
   virtual void apply_A_l_comp_dom_B( TimeEval & time_eval, Cluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double>  & y_out);

};

#endif /* SOLVER_SPECIFIC_ACC_ITERSOLVERACC_H_ */
