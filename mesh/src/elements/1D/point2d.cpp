#include "point2d.h"

void Point2D::normalize()
{
	double l = 1.0 / sqrt(x * x + y * y);
	x *= l;
	y *= l;
}

double Point2D::length() const
{
	return sqrt(x * x + y * y);
}

const Point2D Point2D::operator-() const
{
	return Point2D(-x, -y);
}

void Point2D::flip()
{
	x = -x;
	y = -y;
}

const Point2D Point2D::operator-(const Point2D &p) const
{
	return Point2D(x - p.x, y - p.y);
}

const Point2D Point2D::operator+(const Point2D &p) const
{
	return Point2D(x + p.x, y + p.y);
}

const Point2D Point2D::operator*(double scalar) const
{
	return Point2D(x * scalar, y * scalar);
}

const Point2D Point2D::operator/(double scalar) const
{
	double xscalar = 1. / scalar;
	return Point2D(x * xscalar, y * xscalar);
}

double Point2D::scalar_product_with(const Point2D & p)
{
	return x * p.x + y * p.y;
}

Point2D & Point2D::operator*=(double scalar)
{
	x *= scalar;
	y *= scalar;
	return *this;
}

Point2D & Point2D::operator+=(const Point2D &p)
{
	x += p.x;
	y += p.y;
	return *this;
}

Point2D & Point2D::operator-=(const Point2D &p)
{
	x -= p.x;
	y -= p.y;
	return *this;
}

Point2D & Point2D::operator=(const Point2D &p)
{
	x = p.x;
	y = p.y;
	return *this;
}

Point2D & Point2D::operator/=(double scalar)
{
	x /= scalar;
	y /= scalar;
	return *this;
}
