
#ifndef SRC_MESH_STRUCTURES_REGION_H_
#define SRC_MESH_STRUCTURES_REGION_H_

#include <vector>
#include <map>

namespace espreso {

class Element;
class Evaluator;
class Coordinates;
enum class Property;
enum class ElementType;

struct Region {
	std::string name;
	ElementType eType;
	std::vector<std::map<Property, std::vector<Evaluator*> > > settings;
	mutable double area;

	std::vector<Element*>& elements() { return *_elements; }
	const std::vector<Element*>& elements() const { return *_elements; }

	Region(ElementType eType);
	Region(ElementType eType, std::vector<Element*> &element);

	~Region();

	void computeArea(const Coordinates &coordinates) const;

protected:
	std::vector<Element*> *_elements;
	bool _destroy;
};

}



#endif /* SRC_MESH_STRUCTURES_REGION_H_ */
