
#ifndef SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS_H_
#define SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS_H_

#include "physics.h"
#include "../../basis/point/point.h"

namespace espreso {

struct StructuralMechanicsConfiguration;

struct StructuralMechanics: public virtual Physics
{
	StructuralMechanics(const StructuralMechanicsConfiguration &configuration);

	virtual std::vector<size_t> solutionsIndicesToStore() const;

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const;
	virtual bool isMatrixTimeDependent(const Step &step) const;
	virtual bool isMatrixTemperatureDependent(const Step &step) const;
	virtual void prepare();
	virtual void preprocessData(const Step &step);

	const std::vector<Property>& edgeDOFs() const
	{
		static std::vector<Property> edgeDOFs = { };
		return edgeDOFs;
	}
	const std::vector<Property>& faceDOFs() const
	{
		static std::vector<Property> faceDOFs = { };
		return faceDOFs;
	}
	const std::vector<Property>& elementDOFs() const
	{
		static std::vector<Property> elementDOFs = { };
		return elementDOFs;
	}

protected:
	enum SolutionIndex: size_t {
		DISPLACEMENT = 0,

		SIZE         = 1
	};

	static size_t offset;

	const StructuralMechanicsConfiguration &_configuration;

	// to handle with non-continuous partition
	std::vector<Point> _clusterCenter;
	std::vector<Point> _clusterNorm;
	std::vector<size_t> _clusterNodes;

	std::vector<Point> _domainCenter;
	std::vector<Point> _domainNorm;
	std::vector<size_t> _domainNodes;
};

}


#endif /* SRC_ASSEMBLER_PHYSICS_STRUCTURALMECHANICS_H_ */
