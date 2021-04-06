
#ifndef SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS3DBASE_KERNEL_H_
#define SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS3DBASE_KERNEL_H_

#include "kernel.h"
#include "mover/elasticity.mover.h"

namespace espreso {

class VectorsDense;
class MatrixDense;
struct Element;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
struct Builder;
class HeatTransfer3DKernel;

struct StructuralMechanics3DBaseKernel: public KernelExecutor<StructuralMechanics3DBaseKernel>
{
public:
	StructuralMechanics3DBaseKernel(StructuralMechanics3DBaseKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);
	StructuralMechanics3DBaseKernel(HeatTransfer3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);
	virtual ~StructuralMechanics3DBaseKernel();

	virtual void processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const = 0;
	virtual void processFace(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const = 0;
	virtual void processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const = 0;
	virtual void processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const = 0;

	virtual void elementSolution(ElasticityElementIterator &iterator) = 0;
	virtual void nodeSolution(ElasticityElementIterator &iterator) = 0;

	StructuralMechanicsGlobalSettings &gsettings;
	ElasticityElementIterator iterator;
	std::vector<ElasticityBoundaryIterator> boundaries;
};
}

#endif /* SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS3DBASE_KERNEL_H_ */
