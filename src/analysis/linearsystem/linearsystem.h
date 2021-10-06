
#ifndef SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_

#include "analysis/composer/elementmapping.h"
#include "config/ecf/physics/physicssolver/loadstep.h"
#include "math2/math2.h"
#include "math2/generalization/matrix_base.h"

namespace espreso {

namespace step { struct Step; }

template <typename Assembler, typename Solver = Assembler>
struct AX_LinearSystem {

	template <typename Type>
	struct System {
		Matrix_Base<Type> *A;
		Vector_Base<Type> *x, *b, *dirichlet;
	};

	virtual void setMapping(Matrix_Base<Assembler> *A) const =0;
	virtual void setMapping(Vector_Base<Assembler> *x) const =0;
	virtual void setDirichletMapping(Vector_Base<Assembler> *x) const =0;

	virtual ~AX_LinearSystem() {}

	virtual void info() const {};

	virtual void set(step::Step &step) =0;
	virtual void update(step::Step &step) =0;
	virtual bool solve(step::Step &step) =0;

	System<Assembler> assembler;
	System<Solver> solver;
};



}

#endif /* SRC_ANALYSIS_LINEARSOLVER_LINEARSOLVER_H_ */
