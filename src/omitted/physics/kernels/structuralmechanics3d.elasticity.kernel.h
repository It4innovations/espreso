
#ifndef SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_
#define SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_

#include "structuralmechanics3dbase.kernel.h"
#include "mover/elasticity.mover.h"

namespace espreso {

class VectorsDense;
class MatrixDense;
struct Element;
struct MaterialBaseConfiguration;
struct MaterialConfiguration;
struct Builder;
class HeatTransfer3DKernel;

struct StructuralMechanics3DKernel: public StructuralMechanics3DBaseKernel
{
    StructuralMechanics3DKernel(StructuralMechanics3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);
    StructuralMechanics3DKernel(HeatTransfer3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration);
    virtual ~StructuralMechanics3DKernel();

    virtual void processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const;
    virtual void processFace(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const;
    virtual void processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const;
    virtual void processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const;

    virtual void elementSolution(ElasticityElementIterator &iterator);
    virtual void nodeSolution(ElasticityElementIterator &iterator);
protected:
    void assembleLinearElasticMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, MatrixDense &K, Point &orientation) const;
    void assembleHyperElasticMaterialMatrix(const MaterialBaseConfiguration *mat, MatrixDense &F, MatrixDense &C, MatrixDense &S) const;

    ElementData *orientation;
};

}


#endif /* SRC_PHYSICS_LINEARSYSTEM_KERNELS_STRUCTURALMECHANICS3D_KERNEL_H_ */
