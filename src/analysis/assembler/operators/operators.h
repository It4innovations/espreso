
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_

#include "analysis/assembler/module/assembler.h"

namespace espreso {

struct HeatTransfer;
struct Acoustic;
struct StructuralMechanics;

struct SteadyState;
struct Harmonic;

template <class Setter> void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value, Setter setter);
template <class Setter> void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value, Setter setter);
template <class Setter> void fromExpression2D(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value, Setter setter);
template <class Setter> void fromExpression3D(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value, Setter setter);

void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementNodesValue &value);
void fromExpression(HeatTransfer &module, ParameterData &parameter, ExternalElementGPsValue &value);
void fromExpression(HeatTransfer &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);
void fromExpression(Acoustic &module, ParameterData &parameter, ExternalElementGPsValue &value);
void fromExpression(Acoustic &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);
void fromExpression(StructuralMechanics &module, ParameterData &parameter, ExternalElementNodesValue &value);
void fromExpression(StructuralMechanics &module, ParameterData &parameter, ExternalElementGPsValue &value);
void fromExpression(StructuralMechanics &module, BoundaryParameterPack &parameter, ExternalBoundaryValue &values);

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
void stiffness(HeatTransfer &module);
void RHS(HeatTransfer &module);
void outputGradient(HeatTransfer &module);
void outputFlux(HeatTransfer &module);
void addFiller(HeatTransfer &module, SteadyState &scheme);


void baseFunction(Acoustic &module);
void elementCoordinates(Acoustic &module);
void elementIntegration(Acoustic &module);
void stiffness(Acoustic &module);
void mass(Acoustic &module);
void boundaryMass(Acoustic &module);
void RHS(Acoustic &module);
void addFiller(Acoustic &module, Harmonic &scheme);


void baseFunction(StructuralMechanics &module);
void elementCoordinates(StructuralMechanics &module);
void elementIntegration(StructuralMechanics &module);
void elasticity(StructuralMechanics &module);
void stiffness(StructuralMechanics &module);
void RHS(StructuralMechanics &module);
void addFiller(StructuralMechanics &module, SteadyState &scheme);

}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_OPERATORS_H_ */
