
#ifndef MESH_STRUCTURES_BOUNDARYCONDITION_H_
#define MESH_STRUCTURES_BOUNDARYCONDITION_H_

#include "../elements/elements.h"

namespace espreso {

enum class ConditionElements {
	NODES,
	ELEMENTS
};

enum class ConditionType {
	UNKNOWN,
	DIRICHLET,
	FORCES
};

enum class DOFType : int {
	DISPLACEMENT_X = 0,
	DISPLACEMENT_Y = 1,
	DISPLACEMENT_Z = 2,
	TEMPERATURE = 3,
	PRESSURE = 4,
	TYPES_SIZE = 5
};

class BoundaryCondition {

public:
	virtual BoundaryCondition* copy() =0;

	virtual void fillDirichlet(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<eslocal> &indices, std::vector<double> &values) =0;
	virtual void fillForces(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<double> &f, size_t subdomain) =0;

	void setValue(DOFType type, const std::string &expression)
	{
		_expression[static_cast<int>(type)] = Expression(expression);
		_prescribed[static_cast<int>(type)] = true;
	}

	void setType(ConditionType type)
	{
		_type = type;
	}

	virtual ~BoundaryCondition() {};

protected:
	BoundaryCondition()
	: _type(ConditionType::UNKNOWN),
	  _expression(static_cast<int>(DOFType::TYPES_SIZE), Expression("0")),
	  _prescribed(static_cast<int>(DOFType::TYPES_SIZE), false) {};

	ConditionType _type;
	std::vector<bool> _prescribed;
	std::vector<Expression> _expression;
};

class ElementCondition: public BoundaryCondition {

public:
	std::vector<eslocal> elements;

	void fillDirichlet(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<eslocal> &indices, std::vector<double> &values);
	void fillForces(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<double> &f, size_t subdomain);

	BoundaryCondition* copy()
	{
		return new ElementCondition(*this);
	}
};

class SurfaceCondition: public BoundaryCondition {

public:
	std::vector<Element*> faces;

	void fillDirichlet(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<eslocal> &indices, std::vector<double> &values);
	void fillForces(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<double> &f, size_t subdomain);

	BoundaryCondition* copy()
	{
		return new SurfaceCondition(*this);
	}
	~SurfaceCondition();
};


class NodeCondition: public BoundaryCondition {

public:
	std::vector<eslocal> nodes;
	void fillDirichlet(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<eslocal> &indices, std::vector<double> &values);
	void fillForces(const std::vector<DOFType> &DOFs, const Mesh &mesh, std::vector<double> &f, size_t subdomain);

	BoundaryCondition* copy()
	{
		return new NodeCondition(*this);
	}
};



}

#endif /* MESH_STRUCTURES_BOUNDARYCONDITION_H_ */
