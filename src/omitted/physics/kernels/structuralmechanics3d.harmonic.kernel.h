
#ifndef SRC_PHYSICS_LINEARSYSTEM_KERNELS_HARMONICBALANCE3D_KERNEL_H_
#define SRC_PHYSICS_LINEARSYSTEM_KERNELS_HARMONICBALANCE3D_KERNEL_H_

#include "structuralmechanics3d.elasticity.kernel.h"

namespace espreso {

class VectorsDense;
class MatrixDense;
struct Element;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
struct Builder;


struct HarmonicBalance3DKernel: public StructuralMechanics3DKernel
{
public:
    HarmonicBalance3DKernel(HarmonicBalance3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);
    HarmonicBalance3DKernel(StructuralMechanics3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);
    virtual ~HarmonicBalance3DKernel();

    virtual void processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const;
    virtual void processFace(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &fillere) const;
    virtual void processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const;
    virtual void processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const;

    virtual void elementSolution(ElasticityElementIterator &iterator);
    virtual void nodeSolution(ElasticityElementIterator &iterator);
};
}

#endif /* SRC_PHYSICS_LINEARSYSTEM_KERNELS_HARMONICBALANCE3D_KERNEL_H_ */
