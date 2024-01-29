
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_SUBKERNELS_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_SUBKERNELS_H_

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
#include "analysis/assembler/general/op.thickness.h"
#include "analysis/assembler/general/op.matrix.mass.h"
#include "analysis/assembler/structuralmechanics/op.angularvelocity.h"
#include "analysis/assembler/structuralmechanics/op.acceleration.h"
#include "analysis/assembler/structuralmechanics/op.displacement.h"
#include "analysis/assembler/structuralmechanics/op.smallstraintensor.h"
#include "analysis/assembler/structuralmechanics/op.sigma.h"
#include "analysis/assembler/structuralmechanics/op.stress.h"
#include "analysis/assembler/structuralmechanics/op.elasticity.h"
#include "analysis/assembler/structuralmechanics/op.plasticity.h"
#include "analysis/assembler/structuralmechanics/op.matrix.elasticity.h"
#include "analysis/assembler/structuralmechanics/op.normalpressure.h"
#include "analysis/assembler/structuralmechanics/op.harmonicforce.h"

namespace espreso {

struct StructuralMechanicsOperators {
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
	Elasticity elasticity;
	Plasticity plasticity;
	MatrixElasticity K;
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

struct StructuralMechanicsBoundaryOperators {
	int code;
	size_t elements, chunks;

	size_t esize;
	double surface;

	Basis basis;
	Thickness thickness;
	Coordinates coordinates;
	Integration integration;

	ExternalExpressionVector displacement;
	NormalPressure normalPressure;
	HarmonicForce harmonicForce;

	DataFiller reRHSfiller, imRHSfiller, reDirichlet, imDirichlet;

	struct {
		std::vector<ExternalEvaluator*> node;
		std::vector<ExternalEvaluator*> gp;
	} expressions;
};

}


#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_SUBKERNELS_H_ */
