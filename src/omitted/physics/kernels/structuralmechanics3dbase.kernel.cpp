
#include "structuralmechanics3dbase.kernel.h"
#include "heattransfer3d.kernel.h"
#include "solverdataprovider/structuralmechanics3d.provider.h"
#include "esinfo/meshinfo.h"

using namespace espreso;

StructuralMechanics3DBaseKernel::StructuralMechanics3DBaseKernel(StructuralMechanics3DBaseKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration)
: KernelExecutor(new StructuralMechanics3DSolverDataProvider(configuration)),
  gsettings(gsettings),
  iterator(previous ? &previous->iterator : NULL, physics, gsettings, configuration, 3)
{
    boundaries.reserve(info::mesh->boundaryRegions.size());
    for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
        boundaries.emplace_back(info::mesh->boundaryRegions[i], iterator, configuration, 3);
    }
}

StructuralMechanics3DBaseKernel::StructuralMechanics3DBaseKernel(HeatTransfer3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration)
: KernelExecutor(new StructuralMechanics3DSolverDataProvider(configuration)),
  gsettings(gsettings),
  iterator(previous ? &previous->iterator : NULL, physics, gsettings, configuration, 3)
{
    boundaries.reserve(info::mesh->boundaryRegions.size());
    for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
        boundaries.emplace_back(info::mesh->boundaryRegions[i], iterator, configuration, 3);
    }
}


StructuralMechanics3DBaseKernel::~StructuralMechanics3DBaseKernel()
{
    iterator.clear();
    for (size_t i = 0; i < info::mesh->boundaryRegions.size(); ++i) {
        boundaries[i].clear();
    }
    delete solverDataProvider;
}
