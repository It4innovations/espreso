
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_

#include "analysis/assembler/operator.h"
#include "analysis/assembler/parameter.h"
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"

#include <memory>

namespace espreso {

void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value);
void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value);
void fromExpression(HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);
void fromExpression(Acoustic &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);
void fromExpression(Acoustic &module, ParameterData &parameter, ExternalElementGPsValue &value);


void evaluateFromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value);
void evaluateFromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value);

void averageEnodesToNodes(const ParameterData &from, NodeData &to);
void copyNodesToEnodes(HeatTransfer &module, const NodeData &from, ParameterData &to);
void copyNodesToBnodes(HeatTransfer &module, const NodeData &from, ParameterData &to, size_t region);

void moveEnodesToGPs(HeatTransfer &module, ParameterData &from, ParameterData &to, int dimension);

void baseFunction(HeatTransfer &module);
void elementCoordinates(HeatTransfer &module);
void elementIntegration(HeatTransfer &module);
void thermalConductivity(HeatTransfer &module);
void heatStiffness(HeatTransfer &module);
void heatRHS(HeatTransfer &module);
void outputGradient(HeatTransfer &module);
void outputFlux(HeatTransfer &module);

void addFiller(HeatTransfer &module, SteadyState &scheme);


void baseFunction(Acoustic &module);
void elementCoordinates(Acoustic &module);
void elementIntegration(Acoustic &module);
void acousticStiffness(Acoustic &module);
void acousticMass(Acoustic &module);
void acousticBoundaryMass(Acoustic &module);
void acousticRHS(Acoustic &module);

void addFiller(Acoustic &module, Harmonic &scheme);


}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_ */
