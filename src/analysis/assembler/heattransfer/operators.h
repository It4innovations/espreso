
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OPERATORS_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/basefunctions.h"
#include "analysis/assembler/general/boundarycondition.h"
#include "analysis/assembler/general/material.h"
#include "analysis/assembler/general/filler.h"
#include "analysis/assembler/general/op.coordinates.h"
#include "analysis/assembler/general/op.expression.h"
#include "analysis/assembler/general/op.integration.h"
#include "analysis/assembler/general/op.temperature.h"
#include "analysis/assembler/general/op.thickness.h"
#include "analysis/assembler/general/op.matrix.apply.h"
#include "analysis/assembler/general/op.matrix.mass.h"
#include "analysis/assembler/heattransfer/op.advection.h"
#include "analysis/assembler/heattransfer/op.conductivity.h"
#include "analysis/assembler/heattransfer/op.heatsource.h"
#include "analysis/assembler/heattransfer/op.flux.h"
#include "analysis/assembler/heattransfer/op.gradient.h"
#include "analysis/assembler/heattransfer/op.matrix.conductivity.h"
#include "analysis/assembler/heattransfer/op.externalheat.h"

namespace espreso {

struct HeatTransferElementOperators {
    int code;
    size_t elements, chunks;

    size_t esize;
    double volume;

    double timeIntegrationConstantK, timeIntegrationConstantM;

    Thickness thickness;
    Material material;
    Coordinates coordinates;
    Temperature temperature;
    Integration integration;
    Conductivity conductivity;
    Advection advection;
    BoundaryCondition heatSource, initTemperature;
    MatrixConductivity K;
    MatrixMass M;
    MatrixApply temperatureResidual;

    TemperatureGradient gradient;
    TemperatureFlux flux;

    DataFiller Kfiller, Mfiller, RHSfiller, nRHSfiller;

    struct {
        std::vector<ExternalEvaluator*> node;
        std::vector<ExternalEvaluator*> gp;
    } expressions;
};

struct HeatTransferBoundaryOperators {
    int code;
    size_t elements, chunks;

    size_t esize;
    double surface;

    Thickness thickness;
    Coordinates coordinates;
    Temperature initialTemperature;
    Integration integration;

    ExternalExpression temperature;
    BoundaryCondition heatFlux, heatFlow, htc, externalTemperature;
    ExternalHeat externalHeat;

    DataFiller RHSfiller, dirichlet;

    struct {
        std::vector<ExternalEvaluator*> node;
        std::vector<ExternalEvaluator*> gp;
    } expressions;
};

struct HeatTransferNodeOperators {
    int code;
    size_t elements, chunks;

    size_t esize;

    Thickness thickness;
    Coordinates coordinates;
    Temperature initialTemperature;

    ExternalExpression temperature;
    ExternalHeat externalHeat;

    DataFiller RHSfiller, dirichlet;

    struct {
        std::vector<ExternalEvaluator*> node;
    } expressions;
};


}


#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OPERATORS_H_ */
