
#ifndef ESCONFIG_H_
#define ESCONFIG_H_

#include <cstdlib>

namespace esconfig {

extern int MPIrank;
extern int MPIsize;

namespace solver {
	extern double epsilon;
	extern size_t maxIterations;
}

}


#endif /* ESCONFIG_H_ */
