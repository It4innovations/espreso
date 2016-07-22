
#include "boundarycondition.h"

using namespace espreso;

SurfaceCondition::~SurfaceCondition()
{
	for (size_t i = 0; i < _faces.size(); i++) {
		delete _faces[i];
	}
}

void ElementCondition::set(const std::string &expression, ConditionType type)
{

}

void SurfaceCondition::set(const std::string &expression, ConditionType type)
{

}

void NodeCondition::set(const std::string &expression, ConditionType type)
{
	_expression = Expression(expression);
	_type = type;

	_DOFs.reserve(_nodes.size());
	for (size_t n = 0; n < _nodes.size(); n++) {
		for (int d = 0; d < 3; d++) {
			_DOFs.push_back(3 * _nodes[n] + d);
		}
	}
}


