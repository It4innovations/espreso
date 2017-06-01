
#ifndef SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION_H_
#define SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION_H_

#include "physics.h"

namespace espreso {

enum class Property;
struct AdvectionDiffusionConfiguration;
struct AdvectionDiffusionConvection;

struct AdvectionDiffusion: public virtual Physics
{
	AdvectionDiffusion(const AdvectionDiffusionConfiguration &configuration);

	virtual std::vector<size_t> solutions() const { return { offset + SolutionIndex::TEMPERATURE }; }

	virtual MatrixType getMatrixType(const Step &step, size_t domain) const;
	virtual void prepareTotalFETI();
	virtual void preprocessData(const Step &step);
	virtual void analyticRegularization(size_t domain);

	const std::vector<Property>& pointDOFs() const
	{
		static std::vector<Property> pointDOFs = { Property::TEMPERATURE };
		return pointDOFs;
	}
	const std::vector<Property>& midPointDOFs() const
	{
		static std::vector<Property> midPointDOFs = { Property::TEMPERATURE };
		return midPointDOFs;
	}
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

	virtual ~AdvectionDiffusion() {}

protected:
	enum SolutionIndex: size_t {
		TEMPERATURE = 0,

		SIZE        = 1
	};

	static size_t offset;

	double computeHTC(
			const AdvectionDiffusionConvection &convection, const Element *e, size_t node, Step step,
			double temp) const;

	void convectionMatParameters(
			const AdvectionDiffusionConvection &convection, const Element *e, size_t node, Step step,
			double temp, double T_EXT,
			double &rho, double &dynamic_viscosity, double &dynamic_viscosity_T, double &heat_capacity, double &thermal_conductivity) const;

	const AdvectionDiffusionConfiguration &_configuration;
};

}



#endif /* SRC_ASSEMBLER_PHYSICS_ADVECTIONDIFFUSION_H_ */
