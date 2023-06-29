
#ifndef SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_SUBKERNELLIST_H_
#define SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_SUBKERNELLIST_H_

#include "analysis/assembler/subkernel/basis.h"
#include "analysis/assembler/subkernel/boundarycondition.h"
#include "analysis/assembler/subkernel/thickness.h"
#include "analysis/assembler/subkernel/material.h"
#include "analysis/assembler/subkernel/coordinates.h"
#include "analysis/assembler/subkernel/temperature.h"
#include "analysis/assembler/subkernel/integration.h"
#include "analysis/assembler/subkernel/expression.h"
#include "analysis/assembler/subkernel/filler.h"
#include "analysis/assembler/subkernel/heattransfer/conductivity.h"
#include "analysis/assembler/subkernel/heattransfer/coordinatesystem.h"
#include "analysis/assembler/subkernel/heattransfer/externalheat.h"
#include "analysis/assembler/subkernel/heattransfer/advection.h"
#include "analysis/assembler/subkernel/heattransfer/matrix.h"
#include "analysis/assembler/subkernel/heattransfer/flux.h"
#include "analysis/assembler/subkernel/heattransfer/gradient.h"

namespace espreso {

struct HeatTransferSubKernelsList {
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
	HeatTransferCoordinateSystem coosystem;
	Advection advection;
	BoundaryCondition heatSource;
	HeatTransferMatrix K, M;

	TemperatureGradient gradient;
	TemperatureFlux flux;

	DataFiller Kfiller, Mfiller, RHSfiller, nRHSfiller;

	std::vector<ExternalEvaluator*> expressions;
};

struct HeatTransferBoundarySubKernelsList {
	int code;
	size_t elements, chunks;

	size_t esize;
	double surface;

	Basis basis;
	Thickness thickness;
	Coordinates coordinates;
	Integration integration;

	ExternalExpression temperature;
	BoundaryCondition heatFlux, heatFlow, htc, externalTemperature;
	ExternalHeat externalHeat;

	DataFiller RHSfiller, dirichlet;

	std::vector<ExternalEvaluator*> expressions;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_KERNEL_HEATTRANSFER_SUBKERNELLIST_H_ */
