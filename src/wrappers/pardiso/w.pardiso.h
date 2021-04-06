
#ifndef SRC_WRAPPERS_PARDISO_W_PARDISO_H_
#define SRC_WRAPPERS_PARDISO_W_PARDISO_H_

#include <cstddef>

namespace espreso {

struct PARDISOParameters {
	void* pt[64];
	esint maxfct, mnum, mtype, phase, iparm[64], msglvl, error, *perm;
	double dparm[64];

	int called;

	PARDISOParameters(): pt{}, called(0)
		{
			for (esint i = 0 ; i < 64; i++) {
				pt[i] = NULL;
				iparm[i] = 0;
				dparm[i] = 0;
			}
			maxfct = 1;
			mnum  = 1;
			mtype = 11; // most general for us
			phase = 0;
			msglvl = 0;
			error = 0;
			perm = NULL;
		}
};
}

#endif /* SRC_WRAPPERS_PARDISO_W_PARDISO_H_ */
