
#ifndef SRC_ESINFO_ENVINFO_H_
#define SRC_ESINFO_ENVINFO_H_

namespace espreso {
namespace info {
namespace env {

	extern int MKL_NUM_THREADS;
	extern int OMP_NUM_THREADS;
	extern int SOLVER_NUM_THREADS;
	extern int PAR_NUM_THREADS;

	void set();
	char* pwd();
}
}
}



#endif /* SRC_ESINFO_ENVINFO_H_ */
