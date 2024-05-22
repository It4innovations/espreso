
#ifndef SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OPERATORS_H_
#define SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OPERATORS_H_

#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/general/basis.h"
#include "analysis/assembler/general/boundarycondition.h"
#include "analysis/assembler/general/material.h"
#include "analysis/assembler/general/filler.h"
#include "analysis/assembler/general/op.coordinates.h"
#include "analysis/assembler/general/op.expression.h"
#include "analysis/assembler/general/op.integration.h"
#include "analysis/assembler/general/op.temperature.h"
#include "analysis/assembler/general/op.thickness.h"
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

    Basis basis;
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

    Basis basis;
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

}


#endif /* SRC_ANALYSIS_ASSEMBLER_HEATTRANSFER_OPERATORS_H_ */
