
#ifndef SOLVER_SPECIFIC_CPU_ITERSOLVERCPU_H_
#define SOLVER_SPECIFIC_CPU_ITERSOLVERCPU_H_

#include "feti/specific/itersolver.h"

namespace espreso {

class IterSolverCPU: public IterSolverBase
{
public:
    IterSolverCPU(FETIConfiguration &configuration): IterSolverBase(configuration) {};

    virtual void apply_A_l_comp_dom_B( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);

    virtual void apply_A_l_comp_dom_B_P( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);

    virtual void apply_A_l_comp_dom_B_P_local( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);

    virtual void apply_A_l_comp_dom_B_P_local_sparse( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<esint> & tmp_in_indices, SEQ_VECTOR<double> & tmp_in_values, SEQ_VECTOR<esint> & tmp_out_indices, SEQ_VECTOR<double> & tmp_out_values);

    virtual void Apply_Prec( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );
};

}



#endif /* SOLVER_SPECIFIC_CPU_ITERSOLVERCPU_H_ */
