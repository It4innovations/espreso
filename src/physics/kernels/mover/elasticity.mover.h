
#ifndef SRC_PHYSICS_KERNELS_MOVER_ELASTICITY_MOVER_H_
#define SRC_PHYSICS_KERNELS_MOVER_ELASTICITY_MOVER_H_

#include "mover.h"
#include "moverparameter.h"

#include "esinfo/ecfinfo.h"

namespace espreso {

struct HeatTransferElementIterator;
class PhysicsConfiguration;
class StructuralMechanicsGlobalSettings;
class StructuralMechanicsLoadStepConfiguration;
class StructuralMechanicsConfiguration;

struct ElasticityElementIterator: public ElementIterator {
	MoverCoordinates coordinates;

	MoverFullParameter<InputExpressionOptionalVectorMap, OutputNodes> displacement;
	MoverFullParameter<InputExpressionOptionalVectorMap, OutputNodes> cos;
	MoverFullParameter<InputExpressionOptionalVectorMap, OutputNodes> sin;
	MoverFullParameter<InputExpressionVectorMap        , OutputNodes> acceleration;
	MoverFullParameter<InputExpressionMap              , OutputElements> thickness;

	MoverFullParameter<InputExpressionMap              , OutputNodes> temperature;

	MoverInputParameter<InputExpressionMap> initialTemperature;
//	MoverInputParameter<InputExpressionMap> temperature;
	MoverInputParameter<InputExpressionVectorMap> angularVelocity;

	MoverOutputParameter<OutputNodes> phase;
	MoverOutputParameter<OutputNodes> displacementAmplitude;
	MoverOutputParameter<OutputNodes> velocity;
	MoverOutputParameter<OutputNodes> velocityAmplitude;
	MoverOutputParameter<OutputNodes> accelerationAmplitude;

	MoverOutputParameter<OutputElements> principalStress;
	MoverOutputParameter<OutputElements> componentStress;
	MoverOutputParameter<OutputElements> vonMisesStress;
	MoverOutputParameter<OutputElements> designVariable;
	MoverOutputParameter<OutputElements> complianceDerivation;
	MoverOutputParameter<OutputElements> compliance;

	double minDesignVariable;
	double penaltyFactor;

	bool harmonic, massStabilization;
	bool largeDisplacement;

	StructuralMechanicsLoadStepConfiguration &configuration;
	mutable RotorDynamicsConfiguration::CorotatingRotorConfiguration *corotating;
	mutable RotorDynamicsConfiguration::RotationAxisConfiguration *rotationAxis;
	mutable RotorDynamicsConfiguration::FixedRotorConfiguration *fixed;

	ElasticityElementIterator(ElasticityElementIterator *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration, int dimension, bool omitTemp=false);
	ElasticityElementIterator(HeatTransferElementIterator *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration, int dimension);
};

struct ElasticityBoundaryIterator: public BoundaryIterator {

	struct RotatingForce {
		MoverInputParameter<InputExpressionVector> axis;
		MoverInputParameter<InputExpression> radius;
		MoverInputParameter<InputExpression> mass;
		MoverInputParameter<InputExpression> phaseAngle;
		MoverInputParameter<InputExpression> location;

		RotatingForce(int dimension);
	};

	MoverCoordinates coordinates;
	MoverReducedInputParameter<OutputElements> thickness;

	MoverInputParameter<InputExpression> normalPressure;
	MoverInputParameter<InputExpressionVector> force;
	MoverInputHarmonicParameter harmonicForce;
	RotatingForce rotatingForce;

	ElasticityBoundaryIterator(BoundaryRegionStore *region, ElasticityElementIterator &iterator, StructuralMechanicsLoadStepConfiguration &configuration, int dimension);

	bool hasSettings();
};

}


#endif /* SRC_PHYSICS_KERNELS_MOVER_ELASTICITY_MOVER_H_ */
