
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

class BoundaryCondition {

public:
	const std::vector<eslocal>& DOFs() const
	{
		return _DOFs;
	}

	std::vector<eslocal>& DOFs()
	{
		return _DOFs;
	}

	virtual void set(const std::string &expression, ConditionType type) = 0;

	double value(const Point3D &p)
	{
		return _expression.evaluate(p.x, p.y, p.z);
	}

	ConditionType type()
	{
		return _type;
	}

	virtual ~BoundaryCondition() {};

protected:
	BoundaryCondition(): _type(ConditionType::UNKNOWN), _expression("0") {};

	ConditionType _type;
	Expression _expression;
	std::vector<eslocal> _DOFs;
};

class ElementCondition: public BoundaryCondition {

public:
	ElementCondition(): BoundaryCondition() {};
	BoundaryCondition* copy()
	{
		return new ElementCondition(*this);
	}

	void set(const std::string &expression, ConditionType type);

	const std::vector<eslocal>& elements() const
	{
		return _elements;
	}

	std::vector<eslocal>& elements()
	{
		return _elements;
	}

protected:
	std::vector<eslocal> _elements;

};

class SurfaceCondition: public BoundaryCondition {

public:
	SurfaceCondition(): BoundaryCondition() {};
	BoundaryCondition* copy()
	{
		return new SurfaceCondition(*this);
	}

	void set(const std::string &expression, ConditionType type);

	const std::vector<Element*>& faces() const
	{
		return _faces;
	}

	std::vector<Element*>& faces()
	{
		return _faces;
	}

	~SurfaceCondition();


protected:
	std::vector<Element*> _faces;
};


class NodeCondition: public BoundaryCondition {

public:
	NodeCondition(): BoundaryCondition() {};
	NodeCondition(double value, ConditionType type): BoundaryCondition() {};
	BoundaryCondition* copy()
	{
		return new NodeCondition(*this);
	}

	void set(const std::string &expression, ConditionType type);

	const std::vector<eslocal>& nodes() const
	{
		return _nodes;
	}

	std::vector<eslocal>& nodes()
	{
		return _nodes;
	}

protected:
	std::vector<eslocal> _nodes;
};



}

#endif /* MESH_STRUCTURES_BOUNDARYCONDITION_H_ */
