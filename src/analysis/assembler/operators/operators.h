
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/module/heattransfer.h"

#include <memory>

namespace espreso {

void fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalValue &value);

void baseFunction(AX_HeatTransfer &module);
void elementCoordinates(AX_HeatTransfer &module);
void elementIntegration(AX_HeatTransfer &module);
void thermalConductivity(AX_HeatTransfer &module);
void heatStiffness(AX_HeatTransfer &module);

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_ */
