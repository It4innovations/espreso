#ifndef SOLVER_SPECIFIC_ACC_ITERSOLVERGPU_H_
#define SOLVER_SPECIFIC_ACC_ITERSOLVERGPU_H_

#include "feti/specific/itersolver.h"

namespace espreso {

class IterSolverGPU: public IterSolverBase
{
public:
    IterSolverGPU(FETIConfiguration &configuration): IterSolverBase(configuration) {};

       virtual void apply_A_l_comp_dom_B( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double>  & y_out);
       
       void apply_A_FETI_SC_forFETI   ( SuperCluster & cluster, SEQ_VECTOR<double> & x_in );
       // TODO GPU - tyto funkce se musi upravit na supercluster pro HTFETI 
       void apply_A_FETI_SC_forHFETI  ( Cluster & cluster );  
       void apply_A_htfeti_SC           ( Cluster & cluster );
       // END - TODO_GPU 

       virtual void apply_A_l_comp_dom_B_P             ( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);
       virtual void apply_A_l_comp_dom_B_P_local       ( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);
       virtual void apply_A_l_comp_dom_B_P_local_sparse( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<esint> & tmp_in_indices, SEQ_VECTOR<double> & tmp_in_values, SEQ_VECTOR<esint> & tmp_out_indices, SEQ_VECTOR<double> & tmp_out_values);

       virtual void Apply_Prec( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );


};

}


#endif /* SOLVER_SPECIFIC_ACC_ITERSOLVERGPU_H_ */
