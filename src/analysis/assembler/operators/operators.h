
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"

#include <memory>

namespace espreso {

void fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value);
void fromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value);
void fromExpression(AX_HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);
void fromExpression(AX_Acoustic &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);
void fromExpression(AX_Acoustic &module, ParameterData &parameter, ExternalElementGPsValue &value);


void evaluateFromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value);
void evaluateFromExpression(AX_HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value);

void averageEnodesToNodes(const ParameterData &from, NodeData &to);
void copyNodesToEnodes(AX_HeatTransfer &module, const NodeData &from, ParameterData &to);
void copyNodesToBnodes(AX_HeatTransfer &module, const NodeData &from, ParameterData &to, size_t region);

void moveEnodesToGPs(AX_HeatTransfer &module, ParameterData &from, ParameterData &to, int dimension);

void baseFunction(AX_HeatTransfer &module);
void elementCoordinates(AX_HeatTransfer &module);
void elementIntegration(AX_HeatTransfer &module);
void thermalConductivity(AX_HeatTransfer &module);
void heatStiffness(AX_HeatTransfer &module);
void heatRHS(AX_HeatTransfer &module);
void outputGradient(AX_HeatTransfer &module);
void outputFlux(AX_HeatTransfer &module);

void addFiller(AX_HeatTransfer &module, AX_SteadyState &scheme);


void baseFunction(AX_Acoustic &module);
void elementCoordinates(AX_Acoustic &module);
void elementIntegration(AX_Acoustic &module);
void acousticStiffness(AX_Acoustic &module);
void acousticMass(AX_Acoustic &module);
void acousticBoundaryMass(AX_Acoustic &module);
void acousticRHS(AX_Acoustic &module);

void addFiller(AX_Acoustic &module, AX_Harmonic &scheme);


}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_ */
