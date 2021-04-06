
#ifndef SRC_BASIS_UTILITIES_INLINE_H_
#define SRC_BASIS_UTILITIES_INLINE_H_

#if defined(__INTEL_COMPILER)
#  define FORCE_INLINE __forceinline
#else
#  define FORCE_INLINE inline
#endif

#if defined(__GNUC__)
#  define ALWAYS_INLINE __attribute__((always_inline)) inline
#else
#  define ALWAYS_INLINE FORCE_INLINE
#endif


#endif /* SRC_BASIS_UTILITIES_INLINE_H_ */
