
#ifndef SRC_MESH_STRUCTURES_REGION_H_
#define SRC_MESH_STRUCTURES_REGION_H_

#include <vector>
#include <map>

#include "../settings/property.h"

namespace espreso {

class Element;
class Evaluator;
class Coordinates;

struct Region {
	std::string name;
	std::vector<std::map<Property, std::vector<Evaluator*> > > settings;
	mutable double area;

	std::vector<Element*>& elements() { return *_elements; }
	const std::vector<Element*>& elements() const { return *_elements; }

	Region(): area(0), _destroy(true) { _elements = new std::vector<Element*>(); }
	Region(std::vector<Element*> &element): area(0), _elements(&element), _destroy(false) { }

	~Region()
	{
		if (_destroy) {
			delete _elements;
		}
	}

	void computeArea(const Coordinates &coordinates) const;

protected:
	std::vector<Element*> *_elements;
	bool _destroy;
};

}



#endif /* SRC_MESH_STRUCTURES_REGION_H_ */
