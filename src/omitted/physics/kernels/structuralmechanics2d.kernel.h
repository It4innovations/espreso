
#ifndef SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_
#define SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_

#include "kernel.h"
#include "mover/elasticity.mover.h"

namespace espreso {

class VectorsDense;
class MatrixDense;
struct MaterialBaseConfiguration;
struct Builder;
class HeatTransfer2DKernel;

struct StructuralMechanics2DKernel: public KernelExecutor<StructuralMechanics2DKernel>
{
	StructuralMechanics2DKernel(StructuralMechanics2DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsetting, StructuralMechanicsLoadStepConfiguration &configuration);
	StructuralMechanics2DKernel(HeatTransfer2DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsetting, StructuralMechanicsLoadStepConfiguration &configuration);
	~StructuralMechanics2DKernel();

	void processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const;
	void processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const;
	void processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const;

	void elementSolution(ElasticityElementIterator &iterator);

	ElasticityElementIterator iterator;
	std::vector<ElasticityBoundaryIterator> boundaries;
protected:
	void assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, MatrixDense &K) const;

};

}



#endif /* SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS2D_KERNEL_H_ */
