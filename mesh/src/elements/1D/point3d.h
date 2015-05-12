#ifndef POINT3D_H_
#define POINT3D_H_

#include <math.h>
#include <iostream>
#include <set>

#include <metis.h>

class Point3D
{

public:
	Point3D(): x(0), y(0), z(0) {};
	Point3D(double x, double y, double z): x(x), y(y), z(z) {};
	Point3D(const Point3D &p): x(p.x), y(p.y), z(p.z) {};

	static size_t size() { return 3; }
	void normalize();
	double length() const;
	const Point3D operator-() const;
	void flip();
	const Point3D operator-(const Point3D &) const;
	const Point3D operator+(const Point3D &) const;
	const Point3D operator*(double) const;
	const Point3D operator/(double) const;
	double scalar_product_with(const Point3D &);
	Point3D &operator*=(double);
	Point3D &operator+=(const Point3D &);
	Point3D &operator-=(const Point3D &);
	Point3D &operator=(const Point3D &);
	Point3D &operator/=(double);

	double x, y, z;
};

inline std::ostream& operator<<(std::ostream& os, const Point3D &p)
{
	os << p.x << " " << p.y << " " << p.z;
	return os;
}


inline std::istream& operator>>(std::istream& is, Point3D &p)
{
	if(!(is >> p.x >> p.y >> p.z)) {
		is.setstate(std::ios::failbit);
	}
	return is;
}

inline void operator<<(double* vector, const Point3D &p)
{
	vector[0] = p.x;
	vector[1] = p.y;
	vector[2] = p.z;
}



#endif /* POINT3D_H_ */
