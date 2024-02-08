
#include "fixedwall.h"

#include "esinfo/ecfinfo.h"
#include "math/math.h"

namespace espreso {

template <typename T>
void FixedWall<T>::set(const step::Step &step, FETI<T> &feti)
{
	printf("FIXED WALL\n");
}

template <typename T>
void FixedWall<T>::update(const step::Step &step, FETI<T> &feti)
{

}

template struct FixedWall<double>;

}
