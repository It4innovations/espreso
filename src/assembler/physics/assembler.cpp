
#include "assembler.h"

using namespace espreso;

double Physics::sumNodeProperty(Property property, eslocal node, size_t step, double defaultValue) const
{
	double result = 0;
	bool set = false;
	for (size_t i = 0; i < _mesh.nodes()[node]->regions().size(); i++) {
		if (_mesh.nodes()[node]->regions()[i]->settings.isSet(property)) {
			for (size_t j = 0; j < _mesh.nodes()[node]->regions()[i]->settings[property].size(); j++) {
				result += _mesh.nodes()[node]->regions()[i]->settings[property][j]->evaluate(node);
			}
			set = true;
		}
	}

	return set ? result : defaultValue;
}

double Physics::getNodeProperty(Property property, eslocal node, size_t step, double defaultValue) const
{
	double result = 0;
	bool set = false;
	for (size_t i = 0; i < _mesh.nodes()[node]->regions().size(); i++) {
		if (_mesh.nodes()[node]->regions()[i]->settings.isSet(property)) {
			result = _mesh.nodes()[node]->regions()[i]->settings[property].back()->evaluate(node);
			set = true;
		}
	}

	return set ? result : defaultValue;
}



