
#ifndef SRC_ANALYSIS_ASSEMBLER_KERNEL_STRUCTURALMECHANICS_SUBKERNELLIST_H_
#define SRC_ANALYSIS_ASSEMBLER_KERNEL_STRUCTURALMECHANICS_SUBKERNELLIST_H_

#include "analysis/assembler/subkernel/basis.h"
#include "analysis/assembler/subkernel/boundarycondition.h"
#include "analysis/assembler/subkernel/thickness.h"
#include "analysis/assembler/subkernel/coordinates.h"
#include "analysis/assembler/subkernel/material.h"
#include "analysis/assembler/subkernel/temperature.h"
#include "analysis/assembler/subkernel/integration.h"
#include "analysis/assembler/subkernel/expression.h"
#include "analysis/assembler/subkernel/filler.h"
#include "analysis/assembler/subkernel/structuralmechanics/coordinatesystem.h"
#include "analysis/assembler/subkernel/structuralmechanics/displacement.h"
#include "analysis/assembler/subkernel/structuralmechanics/elasticity.h"
#include "analysis/assembler/subkernel/structuralmechanics/plasticity.h"
#include "analysis/assembler/subkernel/structuralmechanics/matrix.h"
#include "analysis/assembler/subkernel/structuralmechanics/stress.h"

namespace espreso {

struct StructuralMechanicsSubKernelsList {
	int code;
	size_t elements, chunks;

	size_t esize;
	double volume;

	Basis basis;
	Thickness thickness;
	Material material;
	Coordinates coordinates;
	Displacement displacement;
	Temperature temperature;
	Integration integration;
	Elasticity elasticity;
	Plasticity plasticity;
	StructuralMechanicsCoordinateSystem coosystem;
	StructuralMechanicsMatrix K;
	BoundaryCondition acceleration, angularVelocity;
	Stress stress;

	DataFiller Kfiller, RHSfiller, nRHSfiller;

	std::vector<ExternalEvaluator*> expressions;
};

struct StructuralMechanicsBoundarySubKernelsList {
	int code;
	size_t elements, chunks;

	size_t esize;
	double surface;

	Basis basis;
	Thickness thickness;
	Coordinates coordinates;
	Integration integration;

	ExternalExpressionVector displacement;
	BoundaryCondition normalPressure;

	DataFiller RHSfiller, dirichlet;

	std::vector<ExternalEvaluator*> expressions;
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_KERNEL_STRUCTURALMECHANICS_SUBKERNELLIST_H_ */
