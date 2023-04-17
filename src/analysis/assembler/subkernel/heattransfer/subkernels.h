
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_SUBKERNELS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_SUBKERNELS_H_

#include "config/ecf/material/coordinatesystem.h"
#include "math/simd/simd.h"

#include "analysis/assembler/operator.h"
#include "analysis/assembler/module/heattransfer.element.h"

namespace espreso {

struct ConductivityKernel: ActionOperator {
	const char* name() const { return "ConductivityKernel"; }

	ConductivityKernel()
	: conductivity(nullptr), direct(true)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void activate(const ThermalConductivityConfiguration *conductivity, bool direct)
	{
		this->conductivity = conductivity;
		this->direct = direct;
		this->isactive = 1;
	}

	const ThermalConductivityConfiguration *conductivity;
	bool direct;
};

struct HeatTransferCoordinateSystemKernel: ActionOperator {
	const char* name() const { return "HeatTransferCoordinateSystemKernel"; }

	const CoordinateSystemConfiguration *cooSystem;
	CoordinateSystemConfiguration::TYPE type;

	HeatTransferCoordinateSystemKernel()
	: cooSystem(nullptr), type(CoordinateSystemConfiguration::TYPE::CARTESIAN)
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void activate(const CoordinateSystemConfiguration &cooSystem, int isconst)
	{
		this->cooSystem = &cooSystem;
		this->type = cooSystem.type;
		this->isconst = isconst;
		isactive = 1;
		switch (this->type) {
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: this->isconst = 0; break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:   this->isconst = 0; break;
		}
	}
};

struct AdvectionKernel: ActionOperator {
	const char* name() const { return "AdvectionKernel"; }

	ECFExpressionVector *expression;
	double *K;

	AdvectionKernel()
	: expression(nullptr), K(nullptr)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void activate(ECFExpressionVector *expression, double *K)
	{
		this->expression = expression;
		this->K = K;
		if (this->expression) {
			isactive = 1;
		}
	}
};

struct HeatTransferMatrixKernel: public ActionOperator {
	const char* name() const { return "HeatTransferMatrixKernel"; }

	double *K;

	HeatTransferMatrixKernel()
	: K(nullptr)
	{
		isconst = false;
		action = Action::ASSEMBLE | Action::REASSEMBLE;
	}

	void activate(double *K)
	{
		isactive = 1;
		this->K = K;
	}
};

struct TemperatureGradientKernel: ActionOperator {
	const char* name() const { return "TemperatureGradientKernel"; }

	double* gradient;

	TemperatureGradientKernel()
	: gradient(nullptr)
	{
		isconst = false;
		action = Action::SOLUTION;
	}

	void activate(size_t interval, NamedData *gradient)
	{
		this->gradient = gradient->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		isactive = 1;
	}
};

struct TemperatureFluxKernel: ActionOperator {
	const char* name() const { return "TemperatureFluxKernel"; }

	double* flux;

	TemperatureFluxKernel()
	: flux(nullptr)
	{
		isconst = false;
		action = Action::SOLUTION;
	}

	void activate(size_t interval, NamedData *flux)
	{
		this->flux = flux->data.data() + info::mesh->dimension * info::mesh->elements->eintervals[interval].begin;
		isactive = 1;
	}
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_SUBKERNELS_H_ */
