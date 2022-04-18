
#ifndef SRC_MATH2_MATH2_HPP_
#define SRC_MATH2_MATH2_HPP_

#include "math.h"

namespace espreso {
namespace math {

template <typename T> void orthonormalize(Matrix_Dense<T> &m)
{
	for (esint r = 0; r < m.nrows; ++r) {
		for (esint rr = 0; rr < r; ++rr) {
			double scale = math::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + r * m.ncols, 1) / math::dot(m.ncols, m.vals + rr * m.ncols, 1, m.vals + rr * m.ncols, 1);
			math::add(m.ncols, m.vals + r * m.ncols, 1, -scale, m.vals + rr * m.ncols, 1);
		}
		math::scale(m.ncols, 1. / math::norm(m.ncols, m.vals + r * m.ncols, 1), m.vals + r * m.ncols, 1);
	}
}

}
}

#endif /* SRC_MATH2_MATH2_HPP_ */
