
#ifndef PERMONCUBE_H_
#define PERMONCUBE_H_

#include "esmesh.h"

class Permoncube
{

public:
	static void hexahedrons8(
			Mesh &mesh,
			Coordinates &coordinates,
			int *subdomains,
			int *elementsInSub);

	static void tetrahedrons4(
			Mesh &mesh,
			Coordinates &coordinates,
			int *subdomains,
			int *elementsInSub);

	static void dirichlet(
			std::map < int, double >  & dirichlet_x,
			std::map < int, double >  & dirichlet_y,
			std::map < int, double >  & dirichlet_z,
			int *subdomains,
			int *elementsInSub);
};



#endif /* PERMONCUBE_H_ */
