
#include "math/math.h"
#include "math2/math2.h"

#ifdef HAVE_MKL
#include "mkl_service.h"
#endif

using namespace espreso;

void MATH::setNumberOfThreads(int numberOfThreads)
{
#ifdef HAVE_MKL
	mkl_set_num_threads(numberOfThreads);
#endif
}
