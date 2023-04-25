
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_CONDUCTIVITY_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_CONDUCTIVITY_H_

#include "subkernel.h"
#include "config/ecf/material/thermalconductivity.h"

namespace espreso {

struct Conductivity: SubKernel {
	const char* name() const { return "ConductivityKernel"; }

	Conductivity()
	: conductivity(nullptr), indirect(false)
	{
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE;
	}

	void activate(const ThermalConductivityConfiguration *conductivity, bool indirect)
	{
		this->conductivity = conductivity;
		this->indirect = indirect;
		this->isconst = !(this->conductivity->needCoordinates() || this->conductivity->needTemperature());
	}

	const ThermalConductivityConfiguration *conductivity;
	bool indirect;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_HEATTRANSFER_CONDUCTIVITY_H_ */
