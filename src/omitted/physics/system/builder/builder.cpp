
#include "builder.h"

#include "physics/system/wsmpsystem.h"
#include "physics/system/superlusystem.h"
#include "physics/system/pardisosystem.h"
#include "physics/system/mklpdsssystem.h"
#include "physics/system/hypresystem.h"
#include "physics/system/fetisystem.h"
#include "math/matrix.csr.distributed.h"
#include "math/matrix.csr.feti.h"

using namespace espreso;

void Builder::reset(Request matrices, AssemblerData &assembler, SolverData &solver)
{
    if (matrices & Builder::Request::KCM) {
        solver.K->fill(0);
    }
    if (matrices & Builder::Request::R) {
        solver.R->fill(0);
    }
    if (matrices & Builder::Request::f) {
        solver.f->fill(0);
    }
    if (matrices & Builder::Request::BC) {
        solver.BC->fill(0);
    }
}

void Builder::reset(Request matrices, WSMPSystem &system)
{
    reset(matrices, system.assemblers[0], system.solvers[0]);
}

void Builder::reset(Request matrices, SuperLUSystem &system)
{
    reset(matrices, system.assemblers[0], system.solvers[0]);
}

void Builder::reset(Request matrices, PARDISOSystem &system)
{
    reset(matrices, system.assemblers[0], system.solvers[0]);
}

void Builder::reset(Request matrices, MKLPDSSSystem &system)
{
    reset(matrices, system.assemblers[0], system.solvers[0]);
}

void Builder::reset(Request matrices, HYPRESystem &system)
{
    reset(matrices, system.assemblers[0], system.solvers[0]);
}

void Builder::reset(Request matrices, FETISystem &system)
{
    reset(matrices, system.assemblers[0], system.solvers[0]);
}
