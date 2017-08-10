
#ifndef SRC_INPUT_GENERATOR_PRIMITIVES_TRIPLE_H_
#define SRC_INPUT_GENERATOR_PRIMITIVES_TRIPLE_H_

#include <iostream>
#include <cmath>

namespace espreso {
namespace input {

template <typename Ttype>
struct Triple {
	Ttype x, y, z;

	Triple(): x(0), y(0), z(0) {};

	template <typename TOther1, typename TOther2, typename TOther3>
	Triple(TOther1 x, TOther2 y, TOther3 z): x(x), y(y), z(z) {};

	template <typename TOther>
	Triple(TOther* p): x(p[0]), y(p[1]), z(p[2]) {};

	template <typename TOther>
	Triple(const Triple<TOther> &p): x(p.x), y(p.y), z(p.z) {};

	template <typename TOther>
	operator Triple<TOther>() { return Triple<TOther>(x, y, z); }

	template <typename TOther>
	Triple<Ttype> operator+(TOther v) const { return Triple<Ttype>(x + v, y + v, z + v); }

	template <typename TOther>
	Triple<Ttype> operator-(TOther v) const { return Triple<Ttype>(x - v, y - v, z - v); }

	template <typename TOther>
	Triple<Ttype> operator*(TOther v) const { return Triple<Ttype>(x * v, y * v, z * v); }

	template <typename TOther>
	Triple<Ttype> operator/(TOther v) const { return Triple<Ttype>(x / v, y / v, z / v); }

	template <typename TOther>
	Triple<Ttype> operator+(const Triple<TOther> &p) const { return Triple<Ttype>(x + p.x, y + p.y, z + p.z); }

	template <typename TOther>
	Triple<Ttype> operator-(const Triple<TOther> &p) const { return Triple<Ttype>(x - p.x, y - p.y, z - p.z); }

	template <typename TOther>
	Triple<Ttype> operator*(const Triple<TOther> &p) const { return Triple<Ttype>(x * p.x, y * p.y, z * p.z); }

	Triple<double> operator/(const Triple<double> &p) const { return Triple<double>(x / p.x, y / p.y, z / p.z); }

	template <typename TOther>
	Triple<Ttype> operator/(const Triple<TOther> &p) const { return Triple<Ttype>(x / p.x, y / p.y, z / p.z); }

	template <typename TOther>
	Triple<Ttype> operator%(const Triple<TOther> &p) const { return Triple<Ttype>(x % p.x, y % p.y, z % p.z); }

	template <typename TOther>
	Triple<Ttype>& operator+=(const Triple<TOther> &p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}

	bool operator<(Ttype value) const { return x < value || y < value || z < value; }

	template <typename TOther>
	Triple<Ttype>& operator=(const Triple<TOther> &p)
	{
		x = (Ttype)p.x;
		y = (Ttype)p.y;
		z = (Ttype)p.z;
		return *this;
	}

	bool operator==(const Triple<Ttype> &p) { return x == p.x && y == p.y && z == p.z; }

	Ttype sum() const { return x + y + z; }
	Ttype mul() const { return x * y * z; }

	Triple<double> round() const { return Triple<double>(std::round(x), std::round(y), std::round(z)); }
	Triple<double> floor() const { return Triple<double>(std::floor(x), std::floor(y), std::floor(z)); }
	Triple<double> ceil() const { return Triple<double>(std::ceil(x), std::ceil(y), std::ceil(z)); }

	Triple<Ttype> toSize() const { return Triple<Ttype>(1, x, x * y); }
	Triple<Ttype> steps() const { return Triple<Ttype>(x > 1 ? x - 1 : 1, y > 1 ? y - 1 : 1, z > 1 ? z - 1 : 1); }
};

template <typename Ttype>
inline std::ostream& operator<<(std::ostream &os, const Triple<Ttype> &p)
{
	os << p.x << " " << p.y << " " << p.z;
	return os;
}

}
}



#endif /* SRC_INPUT_GENERATOR_PRIMITIVES_TRIPLE_H_ */
