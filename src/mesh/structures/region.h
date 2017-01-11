
#ifndef SRC_MESH_STRUCTURES_REGION_H_
#define SRC_MESH_STRUCTURES_REGION_H_

#include <vector>
#include <map>
#include <string.h>

#include "../settings/property.h"

namespace espreso {

class Element;
class Evaluator;
class Coordinates;

struct Region {
	std::string name;
	std::vector<Element*> elements;
	std::vector<std::map<Property, std::vector<Evaluator*> > > settings;
	mutable double area;

	void computeArea(const Coordinates &coordinates) const;
};

}



#endif /* SRC_MESH_STRUCTURES_REGION_H_ */
