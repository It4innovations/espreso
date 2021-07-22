
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"

#include <memory>

namespace espreso {

void fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementValue &value);
void fromExpression(AX_HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);

void baseFunction(AX_HeatTransfer &module);
void elementCoordinates(AX_HeatTransfer &module);
void elementIntegration(AX_HeatTransfer &module);
void thermalConductivity(AX_HeatTransfer &module);
void heatStiffness(AX_HeatTransfer &module);
void heatRHS(AX_HeatTransfer &module);
void addFiller(AX_HeatTransfer &module);

void baseFunction(AX_Acoustic &module);
void elementCoordinates(AX_Acoustic &module);
void elementIntegration(AX_Acoustic &module);
void acousticStiffness(AX_Acoustic &module);
void acousticMass(AX_Acoustic &module);
void acousticRHS(AX_Acoustic &module);
void addFiller(AX_Acoustic &module);

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_ */
