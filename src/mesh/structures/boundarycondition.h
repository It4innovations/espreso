
#ifndef MESH_STRUCTURES_BOUNDARYCONDITION_H_
#define MESH_STRUCTURES_BOUNDARYCONDITION_H_

#include "../elements/element.h"

namespace espreso {

enum class ConditionElements {
	NODES,
	ELEMENTS
};

enum class ConditionType {
	DIRICHLET,
	FORCES
};

class BoundaryCondition {

public:
	virtual BoundaryCondition* copy() =0;

	const std::vector<eslocal>& DOFs() const
	{
		return _DOFs;
	}

	std::vector<eslocal>& DOFs()
	{
		return _DOFs;
	}

	eslocal DOFsPerNode()
	{
		return _DOFsPerNode;
	}

	double value()
	{
		return _value;
	}

	virtual void set(double value, ConditionType type, const std::vector<bool> &DOFs) = 0;

	ConditionType type()
	{
		return _type;
	}

	virtual ~BoundaryCondition() {};

protected:
	BoundaryCondition(double value, ConditionType type, eslocal DOFsPerNode): _value(value), _type(type), _DOFsPerNode(DOFsPerNode) {};

	std::vector<eslocal> _DOFs;
	eslocal _DOFsPerNode;
	double _value;
	ConditionType _type;
};

class ElementCondition: public BoundaryCondition {

public:
	ElementCondition(): BoundaryCondition(0, ConditionType::DIRICHLET, 3) {};
	BoundaryCondition* copy()
	{
		return new ElementCondition(*this);
	}

	void set(double value, ConditionType type, const std::vector<bool> &DOFs);

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
	SurfaceCondition(): BoundaryCondition(0, ConditionType::DIRICHLET, 3) {};
	BoundaryCondition* copy()
	{
		return new SurfaceCondition(*this);
	}

	void set(double value, ConditionType type, const std::vector<bool> &DOFs);

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
	NodeCondition(): BoundaryCondition(0, ConditionType::DIRICHLET, 3) {};
	NodeCondition(double value, ConditionType type): BoundaryCondition(value, type, 3) {};
	BoundaryCondition* copy()
	{
		return new NodeCondition(*this);
	}

	void set(double value, ConditionType type, const std::vector<bool> &DOFs);

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
