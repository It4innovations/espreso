
#ifndef SRC_BASIS_CONTAINERS__POINT_H_
#define SRC_BASIS_CONTAINERS__POINT_H_

#include <cmath>

namespace espreso {

template <typename TType>
class _Point {

public:
	_Point(): x(0), y(0), z(0) { };
	_Point(TType value): x(value), y(value), z(value) { };
	_Point(TType x, TType y, TType z): x(x), y(y), z(z) { };
	template<typename TOther>
	_Point(const _Point<TOther> &p): x(p.x), y(p.y), z(p.z) { };

	TType& operator[](int i)
	{
		return i == 0 ? x : i == 1 ? y : z;
	}

	const TType& operator[](int i) const
	{
		return i == 0 ? x : i == 1 ? y : z;
	}

	_Point& normalize()
	{
		TType l = 1.0 / length();
		x *= l;
		y *= l;
		z *= l;
		return *this;
	}

	_Point normalize() const
	{
		TType l = 1.0 / length();
		return _Point(x * l, y * l, z * l);
	}

	TType norm() const { return length(); }
	TType length() const { return std::sqrt(x * x + y * y + z * z); }

	const _Point operator-() const
	{
		return _Point(-x, -y, -z);
	}

	void flip()
	{
		x = -x;
		y = -y;
		z = -z;
	}

	void abs()
	{
		x = std::fabs(x);
		y = std::fabs(y);
		z = std::fabs(z);
	}

	_Point operator-(const _Point &p) const
	{
		return _Point(x - p.x, y - p.y, z - p.z);
	}

	_Point operator+(const _Point &p) const
	{
		return _Point(x + p.x, y + p.y, z + p.z);
	}

	_Point operator*(TType scalar) const
	{
		return _Point(x * scalar, y * scalar, z * scalar);
	}

	_Point operator/(TType scalar) const
	{
		TType xscalar = 1. / scalar;
		return _Point(x * xscalar, y * xscalar, z * xscalar);
	}

	TType operator*(const _Point &p) const
	{
		return x * p.x + y * p.y + z * p.z;
	}

	_Point operator/(const _Point &p) const
	{
		return _Point(x / p.x, y / p.y, z / p.z);
	}

	static _Point cross(const _Point &p1, const _Point &p2)
	{
		return _Point(
				p1.y * p2.z - p1.z * p2.y,
				p1.z * p2.x - p1.x * p2.z,
				p1.x * p2.y - p1.y * p2.x);
	}

	static TType cross2d(const _Point &p1, const _Point &p2)
	{
		return p1.x * p2.y - p1.y * p2.x;
	}

	void rodrigues(const _Point &axis, const TType &cos, const TType &sin)
	{
		*this = (*this) * cos + cross(axis, (*this)) * sin + (axis * (axis * (*this))) * (1 - cos);
	}

	_Point &operator*=(TType scalar)
	{
		x *= scalar;
		y *= scalar;
		z *= scalar;
		return *this;
	}

	_Point &operator+=(const _Point &p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}

	_Point &operator-=(const _Point &p)
	{
		x -= p.x;
		y -= p.y;
		z -= p.z;
		return *this;
	}

	_Point &operator/=(TType scalar)
	{
		x /= scalar;
		y /= scalar;
		z /= scalar;
		return *this;
	}

	void getBarycentric(const _Point& p, const _Point& u, const _Point& v, double &s, double &t) const
	{
		_Point w = *this - p;
		double uu = u * u, vv = v * v, uv = u * v, wv = w * v, wu = w * u;
		double d = uv * uv - uu * vv;
		s = (uv * wv - vv * wu) / d;
		t = (uv * wu - uu * wv) / d;
	}

	void minmax(_Point &min, _Point &max)
	{
		min.x = std::fmin(x, min.x);
		min.y = std::fmin(y, min.y);
		min.z = std::fmin(z, min.z);
		max.x = std::fmax(x, max.x);
		max.y = std::fmax(y, max.y);
		max.z = std::fmax(z, max.z);
	}

	void minmax(double *min, double *max)
	{
		min[0] = std::fmin(x, min[0]);
		min[1] = std::fmin(y, min[1]);
		min[2] = std::fmin(z, min[2]);
		max[0] = std::fmax(x, max[0]);
		max[1] = std::fmax(y, max[1]);
		max[2] = std::fmax(z, max[2]);
	}

	TType x, y, z;
};

typedef _Point<double> Point;

template <typename TType>
inline bool operator<(const _Point<TType> &p1, const _Point<TType> &p2) {
	if (p1.x == p2.x) {
		if (p1.y == p2.y) {
			if (p1.z == p2.z) {
				return false;
			}
			return p1.z < p2.z;
		}
		return p1.y < p2.y;
	}
	return p1.x < p2.x;
}

template <typename TType>
inline bool operator>(const _Point<TType> &p1, const _Point<TType> &p2) {
	return p2 < p1;
}

template <typename TType>
inline bool operator==(const _Point<TType> &p1, const _Point<TType> &p2) {
	return !((p1 < p2) || (p2 < p1));
}

template <typename TType>
inline bool operator!=(const _Point<TType> &p1, const _Point<TType> &p2) {
	return ((p1 < p2) || (p2 < p1));
}

}

#endif /* SRC_BASIS_CONTAINERS__POINT_H_ */
