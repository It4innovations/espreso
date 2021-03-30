
#include "structuralmechanics3d.harmonic.kernel.h"

#include "esinfo/stepinfo.h"
#include "mesh/store/nodestore.h"

#include <cmath>

using namespace espreso;


HarmonicBalance3DKernel::HarmonicBalance3DKernel(HarmonicBalance3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanics3DKernel(previous, physics, gsettings, configuration)
{
	if (previous == NULL) {
		solutions.clear();
		solutions.reserve(2);
		solutions.push_back(VectorDense(iterator.cos.output.data->data.size(), iterator.cos.output.data->data.data()));
		solutions.push_back(VectorDense(iterator.sin.output.data->data.size(), iterator.sin.output.data->data.data()));
	} else {
		solutions = previous->solutions;
	}
}

HarmonicBalance3DKernel::HarmonicBalance3DKernel(StructuralMechanics3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanics3DKernel(previous, physics, gsettings, configuration)
{
	solutions.clear();
	solutions.reserve(2);
	solutions.push_back(VectorDense(iterator.cos.output.data->data.size(), iterator.cos.output.data->data.data()));
	solutions.push_back(VectorDense(iterator.sin.output.data->data.size(), iterator.sin.output.data->data.data()));
}

HarmonicBalance3DKernel::~HarmonicBalance3DKernel()
{

}

void HarmonicBalance3DKernel::processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const
{
	StructuralMechanics3DKernel::processElement(builder, iterator, filler);
}

void HarmonicBalance3DKernel::processFace(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	StructuralMechanics3DKernel::processFace(builder, iterator, filler);
}
void HarmonicBalance3DKernel::processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	StructuralMechanics3DKernel::processEdge(builder, iterator, filler);
}
void HarmonicBalance3DKernel::processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	StructuralMechanics3DKernel::processNode(builder, iterator, filler);
}

void HarmonicBalance3DKernel::elementSolution(ElasticityElementIterator &iterator)
{
	if (step::step.type == step::TYPE::FTT) {
		StructuralMechanics3DKernel::elementSolution(iterator);
	}
}
void HarmonicBalance3DKernel::nodeSolution(ElasticityElementIterator &iterator)
{
	if (step::step.type == step::TYPE::FREQUENCY) {
		for (esint n = 0; n < 3; n++) {
			iterator.displacementAmplitude.data[n] = std::sqrt(iterator.cos.data[n] * iterator.cos.data[n] + iterator.sin.data[n] * iterator.sin.data[n]);
			iterator.phase.data[n] = 180 * std::atan2(iterator.sin.data[n], iterator.cos.data[n]) / M_PI;
			iterator.velocityAmplitude.data[n] = step::frequency.angular * iterator.displacementAmplitude.data[n];
			iterator.accelerationAmplitude.data[n] = step::frequency.angular * step::frequency.angular * iterator.displacementAmplitude.data[n];
		}
	}
	if (step::step.type == step::TYPE::FTT) {
		for (esint n = 0; n < 3; n++) {
			double cos = std::cos(step::frequency.angular * step::ftt.time);
			double sin = std::sin(step::frequency.angular * step::ftt.time);
			iterator.displacement.data[n] = iterator.cos.data[n] * cos + iterator.sin.data[n] * sin;
			iterator.velocity.data[n] = -iterator.cos.data[n] * step::frequency.angular * sin + iterator.sin.data[n] * step::frequency.angular * cos;
			iterator.acceleration.data[n] =
					-iterator.cos.data[n] * step::frequency.angular * step::frequency.angular * cos +
					 iterator.sin.data[n] * step::frequency.angular * step::frequency.angular * sin;
		}
	}
}
