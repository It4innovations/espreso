

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include "../../include/espreso.h"

#ifndef esint

#if ESPRESO_LOCAL_INDICES_WIDTH == 64
	typedef long esint;
#elif ESPRESO_LOCAL_INDICES_WIDTH == 32
	typedef int esint;
#else
	#error "Incorrect user-supplied value for ESPRESO_LOCAL_INDICES_WIDTH"
#endif

#endif

#ifndef eslong

#if ESPRESO_GLOBAL_INDICES_WIDTH == 64
	typedef long eslong;
#elif ESPRESO_GLOBAL_INDICES_WIDTH == 32
	typedef int eslong;
#else
	#error "Incorrect user-supplied value for ESPRESO_GLOBAL_INDICES_WIDTH"
#endif

#endif


#endif /* DEFINITIONS_H_ */
