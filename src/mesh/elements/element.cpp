
#include <fstream>

#include "element.h"
#include "../structures/region.h"
#include "../structures/coordinates.h"
#include "../settings/evaluator.h"

using namespace espreso;

void Element::store(std::ofstream& os, const Coordinates &coordinates, size_t part)
{
	eslocal value = vtkCode(), pSize = params(), p;
	os.write(reinterpret_cast<const char *>(&value), sizeof(eslocal));
	for (size_t n = 0; n < nodes(); n++) {
		eslocal index = coordinates.localIndex(node(n), part);
		os.write(reinterpret_cast<const char *>(&index), sizeof(eslocal));
	}
	os.write(reinterpret_cast<const char *>(&pSize), sizeof(eslocal));
	if (pSize) {
		p = param(Element::MATERIAL);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
		p = param(Element::CONSTANT);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
		p = param(Element::COORDINATES);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
		p = param(Element::BODY);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
	}
}

void Element::rotateOutside(const Element* parent, const Coordinates &coordinates, Point &normal) const
{
	Point eMid(0, 0, 0), mid(0, 0, 0);
	for (size_t i = 0; i < parent->coarseNodes(); i++) {
		eMid += coordinates[parent->node(i)];
	}
	eMid /= parent->coarseNodes();

	for (size_t i = 0; i < coarseNodes(); i++) {
		mid += coordinates[node(i)];
	}
	mid /= coarseNodes();

	Point outside = mid - eMid;
	if (outside.x * normal.x + outside.y * normal.y + outside.z * normal.z > 0) {
		normal.flip();
	}
}

bool Element::hasProperty(Property property, size_t step) const
{
	for (size_t i = 0; i < _regions.size(); i++) {
		if (step < _regions[i]->settings.size() && _regions[i]->settings[step].count(property)) {
			return true;
		}
	}
	return false;
}

double Element::sumProperty(Property property, eslocal node, size_t step, double defaultValue) const
{
	double result = 0;
	bool set = false;
	for (size_t i = 0; i < _regions.size(); i++) {
		if (step < _regions[i]->settings.size()) {
			auto it = _regions[i]->settings[step].find(property);
			if (it != _regions[i]->settings[step].end()) {
				for (size_t j = 0; j < it->second.size(); j++) {
					result += it->second[j]->evaluate(this->node(node));
				}
				set = true;
			}
		}
	}

	return set ? result : defaultValue;
}

double Element::getProperty(Property property, eslocal node, size_t step, double defaultValue) const
{
	for (size_t i = 0; i < _regions.size(); i++) {
		if (step < _regions[i]->settings.size()) {
			auto it = _regions[i]->settings[step].find(property);
			if (it == _regions[i]->settings[step].end()) {
				continue;
			}
			return it->second.back()->evaluate(this->node(node));
		}
	}

	return defaultValue;
}



