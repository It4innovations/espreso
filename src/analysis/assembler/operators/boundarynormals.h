#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_BOUNDARYNORMALS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_BOUNDARYNORMALS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/math.hpp"

#include <iostream>

namespace espreso {

template <size_t nodes, size_t gps>
struct BoundaryNormal2D: public ActionOperator {

	InputParameterIterator dN, coords;
	OutputParameterIterator normal;

	BoundaryNormal2D(
		int interval,
		const ParameterData &dN,
		const ParameterData &coords,
		ParameterData &normal)
	: dN(dN, interval),
	  coords(coords, interval),
	  normal(normal, interval)
	{

	}

	void operator++()
	{
		++dN;
		++coords;
		++normal;
	}

	void move(int n)
	{
		dN += n;
		coords += n;
		normal += n;
	}

	void operator()()
	{
		double tangent[2];
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			M1NMN2<nodes>(1.0, dN.data + nodes * gpindex, coords.data, tangent);
			double norm = std::sqrt(tangent[0]*tangent[0] + tangent[1]*tangent[1]);
			tangent[0] /= norm;
			tangent[1] /= norm;
			normal[2*gpindex + 0] = tangent[1];
			normal[2*gpindex + 1] = -tangent[0];
		}
	}
};

template <size_t nodes, size_t gps>
struct BoundaryNormal3D: public ActionOperator {

	InputParameterIterator dN, coords;
	OutputParameterIterator normals;

	BoundaryNormal3D(
		int interval,
		const ParameterData &dN,
		const ParameterData &coords,
		ParameterData &normals)
	: dN(dN, interval),
	  coords(coords, interval),
	  normals(normals, interval)
	{

	}

	void operator++()
	{
		++dN;
		++coords;
		++normals;
	}

	void move(int n)
	{
		dN += n;
		coords += n;
		normals += n;
	}

	void operator()()
	{
		double tangents[6];
		for (size_t gpindex = 0; gpindex < gps; ++gpindex) {
			M2NMN3<nodes>(1.0, dN.data + 2 * nodes * gpindex, coords.data, tangents);
			
			Point v1(tangents[0], tangents[1], tangents[2]);
			Point v2(tangents[3], tangents[4], tangents[5]);
			Point normal = Point::cross(v1, v2);
			normal.normalize();

			// std::cout << "    >>>   " << normal.x << "  " << normal.y << "  " << normal.z << "\n";
			
			normals[3*gpindex + 0] = normal[0];
			normals[3*gpindex + 1] = normal[1];
			normals[3*gpindex + 2] = normal[2];
 		}
	}
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_BOUNDARYNORMALS_H_ */
