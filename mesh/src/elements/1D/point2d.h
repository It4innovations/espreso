#ifndef POINT2D_H_
#define POINT2D_H_

#include <math.h>
#include <iostream>
#include <set>

#include <metis.h>

class Point2D
{

public:
	Point2D(): x(0), y(0) {};
	Point2D(real_t x, real_t y): x(x), y(y) {};
	Point2D(const Point2D &p): x(p.x), y(p.y) {};

	static size_t size() { return 2; }
	void normalize();
	real_t length() const;
	const Point2D operator-() const;
	void flip();
	const Point2D operator-(const Point2D &) const;
	const Point2D operator+(const Point2D &) const;
	const Point2D operator*(real_t) const;
	const Point2D operator/(real_t) const;
	real_t scalar_product_with(const Point2D &);
	Point2D &operator*=(real_t);
	Point2D &operator+=(const Point2D &);
	Point2D &operator-=(const Point2D &);
	Point2D &operator=(const Point2D &);
	Point2D &operator/=(real_t);

	real_t x, y;
};

inline std::ostream& operator<<(std::ostream& os, const Point2D &p)
{
	os << p.x << " " << p.y;
	return os;
}

inline std::istream& operator>>(std::istream& is, Point2D &p)
{
	if(!(is >> p.x >> p.y)) {
		is.setstate(std::ios::failbit);
	}
	return is;
}

inline void operator<<(double* vector, const Point2D &p)
{
	vector[0] = p.x;
	vector[1] = p.y;
}

#endif /* POINT2D_H_ */
