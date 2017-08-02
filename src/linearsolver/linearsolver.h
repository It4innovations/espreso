
#ifndef SRC_LINEARSOLVER_LINEARSOLVER_H_
#define SRC_LINEARSOLVER_LINEARSOLVER_H_

namespace espreso {

enum Matrices: int;

class LinearSolver {

public:
	virtual void update(Matrices matrices) =0;
	virtual void solve() =0;
	virtual void finalize() =0;

	virtual ~LinearSolver() {}
};

}



#endif /* SRC_LINEARSOLVER_LINEARSOLVER_H_ */
