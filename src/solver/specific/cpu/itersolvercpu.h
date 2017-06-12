
#ifndef SOLVER_SPECIFIC_CPU_ITERSOLVERCPU_H_
#define SOLVER_SPECIFIC_CPU_ITERSOLVERCPU_H_

#include "../itersolver.h"

namespace espreso {

class IterSolverCPU: public IterSolverBase
{
public:
	IterSolverCPU(const ESPRESOSolver &configuration): IterSolverBase(configuration) {};

    virtual void apply_A_l_comp_dom_B( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out);

    virtual void Apply_Prec( TimeEval & time_eval, SuperCluster & cluster, SEQ_VECTOR<double> & x_in, SEQ_VECTOR<double> & y_out );
};

}



#endif /* SOLVER_SPECIFIC_CPU_ITERSOLVERCPU_H_ */
