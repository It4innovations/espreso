
#include "esinfo/eslog.h"
#include "math/wrappers/math.lapack.h"

#include <complex>

#ifndef HAVE_MKL
#ifdef HAVE_LAPACK
#include "lapacke.h"

namespace espreso {
namespace math {
namespace lapack {

}
}
}

#endif
#endif
