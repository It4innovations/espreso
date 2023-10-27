
#include "simd.sve.h"

#if defined(__ARM_FEATURE_SVE)
__svep SIMD::mask = svptrue_b64();
#endif
