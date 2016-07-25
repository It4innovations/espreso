
#ifndef SRC_MESH_STRUCTURES_INITIALCONDITION_H_
#define SRC_MESH_STRUCTURES_INITIALCONDITION_H_

#include "../elements/element.h"
#include <string>
#include <map>

namespace espreso {

class InitialCondition {

public:
	virtual InitialCondition* copy() =0;

	void setValue(const std::string &name, const std::string &expression)
	{
		_expression[name] = new Expression(expression);
	}

	virtual ~InitialCondition()
	{
		for (auto it = _expression.begin(); it != _expression.end(); ++it) {
			delete it->second;
		}
	};

protected:
	std::map<std::string, Expression*> _expression;
};

class ElementInitialization: public InitialCondition {

public:
	std::vector<eslocal> elements;


	ElementInitialization* copy()
	{
		return new ElementInitialization(*this);
	}
};

class SurfaceInitialization: public InitialCondition {

public:
	std::vector<Element*> faces;


	SurfaceInitialization* copy()
	{
		return new SurfaceInitialization(*this);
	}
	~SurfaceInitialization();
};


class NodeInitialization: public InitialCondition {

public:
	std::vector<eslocal> nodes;

	NodeInitialization* copy()
	{
		return new NodeInitialization(*this);
	}
};

}



#endif /* SRC_MESH_STRUCTURES_INITIALCONDITION_H_ */
