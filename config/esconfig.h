
#ifndef ESCONFIG_H_
#define ESCONFIG_H_

#include <cstdlib>

namespace esconfig {

enum Discretization {
	FEM,
	BEM
};

extern int MPIrank;
extern int MPIsize;
extern Discretization discretization;

namespace solver {
	extern double epsilon;
	extern size_t maxIterations;
}

}


#endif /* ESCONFIG_H_ */
