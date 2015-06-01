
#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include "../../include/espreso.h"

#if ESPRESO_LOCAL_INDICES_WIDTH == 64
	typedef long long esint;
#else
	typedef int esint;
#endif

#if ESPRESO_GLOBAL_INDICES_WIDTH == 32
	typedef int eslong;
#else
	typedef long long eslong;
#endif



#endif /* DEFINITIONS_H_ */
