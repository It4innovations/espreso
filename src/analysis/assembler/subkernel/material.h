
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_MATERIAL_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_MATERIAL_H_

#include "subkernel.h"
#include "config/ecf/material/material.h"

namespace espreso {

struct Material: SubKernel {
	const char* name() const { return "Material"; }

	Material()
	: configuration(nullptr)
	{
		action = Assembler::ASSEMBLE | Assembler::REASSEMBLE;
	}

	void activate(const MaterialConfiguration *configuration)
	{
		this->configuration = configuration;
	}

	const MaterialConfiguration *configuration;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_MATERIAL_H_ */
