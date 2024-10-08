
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OPERATORS_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/basis.h"
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
#include "analysis/assembler/structuralmechanics/op.angularvelocity.h"
#include "analysis/assembler/structuralmechanics/op.acceleration.h"
#include "analysis/assembler/structuralmechanics/op.displacement.h"
#include "analysis/assembler/structuralmechanics/op.smallstraintensor.h"
#include "analysis/assembler/structuralmechanics/op.sigma.h"
#include "analysis/assembler/structuralmechanics/op.stress.h"
#include "analysis/assembler/structuralmechanics/op.linearelasticity.h"
#include "analysis/assembler/structuralmechanics/op.hyperelasticity.h"
#include "analysis/assembler/structuralmechanics/op.plasticity.h"
#include "analysis/assembler/structuralmechanics/op.matrix.linearelasticity.h"
#include "analysis/assembler/structuralmechanics/op.matrix.hyperelasticity.h"
#include "analysis/assembler/structuralmechanics/op.matrix.largedisplacement.h"
#include "analysis/assembler/structuralmechanics/op.matrix.corotation.h"
#include "analysis/assembler/structuralmechanics/op.normalpressure.h"
#include "analysis/assembler/structuralmechanics/op.pressure.h"
#include "analysis/assembler/structuralmechanics/op.harmonicforce.h"
#include "analysis/assembler/structuralmechanics/op.fluidforce.h"

namespace espreso {

struct StructuralMechanicsElementOperators {
    int code;
    size_t elements, chunks;

    size_t esize;
    double volume;

    Basis basis;
    Thickness thickness;
    Material material;
    Coordinates coordinates;
    Displacement displacement;
    SmallStrainTensor smallStrainTensor;
    Temperature temperature;
    Integration integration;
    LinearElasticity linearElasticity;
    HyperElasticity hyperElasticity;
    Plasticity plasticity;
    MatrixLinearElasticity matrixLinearElasticity;
    MatrixHyperElasticity matrixHyperElasticity;
    MatrixLargeDisplacement largeDisplacement;
    MatrixCorotation corotation;
    MatrixMass M;
    ElementCondition acceleration, angularVelocity;
    Sigma sigma;
    Stress stress;

    DataFiller Kfiller, Mfiller, Cfiller, reRHSfiller, imRHSfiller, reNRHSfiller, imNRHSfiller;

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

    Basis basis;
    Thickness thickness;
    Coordinates coordinates;
    Integration integration;
    Displacement displacement;
    Normal normal;

    Pressure pressure;
    NormalPressure normalPressure;
    FluidForce fluidForce;

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

    Basis basis;
    Thickness thickness;
    Coordinates coordinates;
    Normal normal;

    ExternalExpressionVector displacement;
    HarmonicForce harmonicForce;

    DataFiller reRHSfiller, imRHSfiller, reDirichlet, imDirichlet;

    struct {
        std::vector<ExternalEvaluator*> node;
    } expressions;
};

}


#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_OPERATORS_H_ */
