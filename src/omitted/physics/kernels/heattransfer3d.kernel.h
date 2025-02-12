
#ifndef SRC_PHYSICS_LINEARSYSTEM_KERNELS_HEATTRANSFER3D_KERNEL_H_
#define SRC_PHYSICS_LINEARSYSTEM_KERNELS_HEATTRANSFER3D_KERNEL_H_

#include "kernel.h"
#include "mover/heattransfer.mover.h"

namespace espreso {

class MatrixDense;
class VectorsDense;
struct MaterialBaseConfiguration;
struct Builder;

struct HeatTransfer3DKernel: public KernelExecutor<HeatTransfer3DKernel>
{
    HeatTransfer3DKernel(HeatTransfer3DKernel *previous, PhysicsConfiguration &physics, HeatTransferGlobalSettings &gsettings, HeatTransferLoadStepConfiguration &configuration);
    ~HeatTransfer3DKernel();

    void processElement(const Builder &builder, const HeatTransferElementIterator &iterator, InstanceFiller &output) const;
    void processFace(const Builder &builder, const HeatTransferBoundaryIterator &iterator, InstanceFiller &output) const;
    void processEdge(const Builder &builder, const HeatTransferBoundaryIterator &iterator, InstanceFiller &output) const;

    void elementSolution(const HeatTransferElementIterator &iterator);

    HeatTransferGlobalSettings &gsettings;
    HeatTransferElementIterator iterator;
    std::vector<HeatTransferBoundaryIterator> boundaries;
protected:
    void assembleMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double phase, double time, double temp, MatrixDense &K, MatrixDense &CD, bool tangentCorrection) const;

};

}


#endif /* SRC_PHYSICS_LINEARSYSTEM_KERNELS_HEATTRANSFER3D_KERNEL_H_ */
