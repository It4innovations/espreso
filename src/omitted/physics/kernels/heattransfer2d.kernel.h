
#ifndef SRC_PHYSICS_LINEARSYSTEM_KERNELS_HEATTRANSFER2D_KERNEL_H_
#define SRC_PHYSICS_LINEARSYSTEM_KERNELS_HEATTRANSFER2D_KERNEL_H_

#include "kernel.h"
#include "mover/heattransfer.mover.h"

namespace espreso {

class VectorsDense;
class MatrixDense;
struct MaterialBaseConfiguration;
struct Builder;

class HeatTransfer2DKernel: public KernelExecutor<HeatTransfer2DKernel>
{
public:
	HeatTransfer2DKernel(HeatTransfer2DKernel *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);
	~HeatTransfer2DKernel();

	void processElement(const Builder &builder, HeatTransferElementIterator &iterator, InstanceFiller &output) const;
	void processEdge(const Builder &builder, HeatTransferBoundaryIterator &iterator, InstanceFiller &output) const;

	void elementSolution(HeatTransferElementIterator &iterator);

	HeatTransferGlobalSettings &gsettings;
	HeatTransferElementIterator iterator;
	std::vector<HeatTransferBoundaryIterator> boundaries;

protected:
	void assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, MatrixDense &K, MatrixDense &CD, bool tangentCorrection) const;

};

}


#endif /* SRC_PHYSICS_LINEARSYSTEM_KERNELS_HEATTRANSFER2D_KERNEL_H_ */
