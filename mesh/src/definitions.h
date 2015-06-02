
#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include "../../include/espreso.h"

#if ESPRESO_LOCAL_INDICES_WIDTH == 64
	typedef long long esint;
#elif ESPRESO_LOCAL_INDICES_WIDTH == 32
	typedef int esint;
#else
	#error "Incorrect user-supplied value for ESPRESO_LOCAL_INDICES_WIDTH"
#endif

#if ESPRESO_GLOBAL_INDICES_WIDTH == 64
	typedef long long eslong;
#elif ESPRESO_GLOBAL_INDICES_WIDTH == 32
	typedef int eslong;
#else
	#error "Incorrect user-supplied value for ESPRESO_GLOBAL_INDICES_WIDTH"
#endif



#endif /* DEFINITIONS_H_ */
