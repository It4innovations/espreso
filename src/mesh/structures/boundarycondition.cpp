
#include "boundarycondition.h"

using namespace espreso;

SurfaceCondition::~SurfaceCondition()
{
	for (size_t i = 0; i < _faces.size(); i++) {
		delete _faces[i];
	}
}

void ElementCondition::set(double value, ConditionType type, const std::vector<bool> &DOFs)
{

}

void SurfaceCondition::set(double value, ConditionType type, const std::vector<bool> &DOFs)
{

}

void NodeCondition::set(double value, ConditionType type, const std::vector<bool> &DOFs)
{
	_value = value;
	_type = type;

	size_t size = 0;
	for (int d = 0; d < DOFs.size(); d++) {
		if (DOFs[d]) {
			size++;
		}
	}
	size *= _nodes.size();

	_DOFs.reserve(size);
	for (size_t n = 0; n < _nodes.size(); n++) {
		for (int d = 0; d < DOFs.size(); d++) {
			if (DOFs[d]) {
				_DOFs.push_back(DOFs.size() * _nodes[n] + d);
			}
		}
	}
}


