
#include "assembler.h"

using namespace espreso;

double Physics::sumNodeProperty(Property property, eslocal node, size_t step, double defaultValue) const
{
	double result = 0;
	bool set = false;
	for (size_t i = 0; i < _mesh.nodes()[node]->regions().size(); i++) {
		if (step < _mesh.nodes()[node]->regions()[i]->settings.size()) {
			auto it = _mesh.nodes()[node]->regions()[i]->settings[step].find(property);
			if (it == _mesh.nodes()[node]->regions()[i]->settings[step].end()) {
				continue;
			}
			for (size_t j = 0; j < it->second.size(); j++) {
				result += it->second[j]->evaluate(node);
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
		if (step < _mesh.nodes()[node]->regions()[i]->settings.size()) {
			auto it = _mesh.nodes()[node]->regions()[i]->settings[step].find(property);
			if (it == _mesh.nodes()[node]->regions()[i]->settings[step].end()) {
				continue;
			}
			result = it->second.back()->evaluate(node);
			set = true;
		}
	}

	return set ? result : defaultValue;
}



