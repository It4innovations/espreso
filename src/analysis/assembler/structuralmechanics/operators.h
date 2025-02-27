
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OPERATORS_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/basefunctions.h"
#include "analysis/assembler/general/boundarycondition.h"
#include "analysis/assembler/general/material.h"
#include "analysis/assembler/general/filler.h"
#include "analysis/assembler/general/op.coordinates.h"
#include "analysis/assembler/general/op.expression.h"
#include "analysis/assembler/general/op.integration.h"
#include "analysis/assembler/general/op.normal.h"
#include "analysis/assembler/general/op.temperature.h"
#include "analysis/assembler/general/op.thickness.h"
#include "analysis/assembler/general/op.thickness.h"
#include "analysis/assembler/general/op.matrix.mass.h"
#include "analysis/assembler/general/op.print.eigenvalues.h"
#include "analysis/assembler/structuralmechanics/op.angularvelocity.h"
#include "analysis/assembler/structuralmechanics/op.acceleration.h"
#include "analysis/assembler/structuralmechanics/op.displacement.h"
#include "analysis/assembler/structuralmechanics/op.stress.h"
#include "analysis/assembler/structuralmechanics/op.material.h"
#include "analysis/assembler/structuralmechanics/op.matrix.elasticity.h"
#include "analysis/assembler/structuralmechanics/op.matrix.corotation.h"
#include "analysis/assembler/structuralmechanics/op.normalpressure.h"
#include "analysis/assembler/structuralmechanics/op.pressure.h"
#include "analysis/assembler/structuralmechanics/op.force.h"
#include "analysis/assembler/structuralmechanics/op.harmonicforce.h"
#include "analysis/assembler/structuralmechanics/op.velocity.h"

namespace espreso {

struct StructuralMechanicsElementOperators {
    int code;
    size_t elements, chunks;

    size_t esize;
    double volume;

    Thickness thickness;
    Coordinates coordinates;
    Displacement displacement;
    Temperature temperature;
    Velocity velocity;
    Integration integration;
    MaterialStructuralMechanics material;
    BoundaryCondition initVelocity;
    MatrixElasticity matrixElasticity;
    MatrixCorotation matrixCorotation;
    MatrixMass M, postM;
    ElementCondition acceleration, angularVelocity;
    Stress stress;
    PrintEigenValues print;

    DataFiller Kfiller, Mfiller, Cfiller, reRHSfiller, imRHSfiller, reNRHSfiller, imNRHSfiller;
    DataFiller postMfiller, postBfiller;

    struct {
        std::vector<ExternalEvaluator*> node;
        std::vector<ExternalEvaluator*> gp;
    } expressions;
};

struct StructuralMechanicsFaceOperators {
    int code;
    size_t elements, chunks;

    size_t esize;
    double surface;

    Thickness thickness;
    Coordinates coordinates;
    Integration integration;
    Displacement displacement;
    Normal normal;

    Force force;
    Pressure pressure;
    NormalPressure normalPressure;

    DataFiller reRHSfiller, imRHSfiller;

    struct {
        std::vector<ExternalEvaluator*> node;
        std::vector<ExternalEvaluator*> gp;
    } expressions;
};

struct StructuralMechanicsNodeOperators {
    int code;
    size_t elements, chunks;

    size_t esize;

    Thickness thickness;
    Coordinates coordinates;
    Normal normal;

    ExternalExpressionVector displacement;
    Force force;
    HarmonicForce harmonicForce;

    DataFiller reRHSfiller, imRHSfiller, reDirichlet, imDirichlet;

    struct {
        std::vector<ExternalEvaluator*> node;
    } expressions;
};

}


#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OPERATORS_H_ */
