
#ifndef SRC_BASIS_SFC_HILBERTCURVE_H_
#define SRC_BASIS_SFC_HILBERTCURVE_H_

#include "spacefillingcurve.h"

namespace espreso {

template <typename T>
struct HilbertCurve: public SpaceFillingCurve<T> {

	HilbertCurve(size_t dimension, size_t depth, size_t npoints, _Point<T>* coordinates): SpaceFillingCurve<T>(dimension, depth, npoints, coordinates) {};

	void D1toD2(size_t n, size_t d, size_t &x, size_t &y) const
	{
		int rx, ry;
		x = y = 0;
		for (size_t s = 1, t = d; s < n; s *= 2, t /= 4) {
			rx = 1 & (t / 2);
			ry = 1 & (t ^ rx);

			rotateD2(s, x, y, rx, ry);

			x += s * rx;
			y += s * ry;
		}
	}

	void D1toD3(size_t n, size_t d, size_t &x, size_t &y, size_t &z) const
	{
		int rx, ry, rz;
		x = y = z = 0;
		for (size_t s = 1, t = d; s < n; s *= 2, t /= 8) {
			rx = (1 & (t / 4)) ^ (1 & (t / 2));
			ry = (1 & (t / 4));
			rz = (1 & (t / 2)) ^ (1 & t);

			rotateD3(s, x, y, z, rx, ry, rz);

			x += s * rx;
			y += s * ry;
			z += s * rz;
		}
	}

	size_t D2toD1(size_t n, size_t x, size_t y) const
	{
		int rx, ry;
		size_t d = 0;

		for (size_t s = n / 2; s > 0; s /= 2) {
			rx = (x & s) > 0;
			ry = (y & s) > 0;
			d += s * s * ((3 * rx) ^ ry);

			rotateD2(s, x, y, rx, ry);
		}
		return d;
	}

	size_t D3toD1(size_t n, size_t x, size_t y, size_t z) const
	{
		int rx, ry, rz;
		size_t d = 0;

		for (size_t s = n / 2; s > 0; s /= 2) {
			rx = (x & s) > 0;
			ry = (y & s) > 0;
			rz = (z & s) > 0;
			d += s * s * s * (4 * ry + ((3 * (ry ^ rx)) ^ rz));

			rotateD3(s, x, y, z, rx, ry, rz);
		}
		return d;
	}

private:
	void rotateD2(size_t n, size_t &x, size_t &y, int rx, int ry) const
	{
		if (ry == 0) {
			if (rx == 1) {
				x = n - 1 - x;
				y = n - 1 - y;
			}
			std::swap(x, y);
		}
	}
	void rotateD3(size_t n, size_t &x, size_t &y, size_t &z, int rx, int ry, int rz) const
	{
		if (rz == 0) {
			if (rx == 1) {
				std::swap(x, z);
				x = n - 1 - x;
				z = n - 1 - z;
			} else {
				if (ry == 1) {
					y = n - 1 - y;
					z = n - 1 - z;
				}
				std::swap(y, z);
			}
		} else {
			std::swap(x, y);
			if (ry == 1) {
				x = n - 1 - x;
				y = n - 1 - y;
			}
		}
	}
};

}

#endif /* SRC_BASIS_SFC_HILBERTCURVE_H_ */
