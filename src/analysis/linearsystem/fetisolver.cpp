
#include "fetisolver.h"

#include "regularization/regularization.h"
#include "constrains/constrains.h"

#include "analysis/physics/heat.steadystate.linear.h"
#include "analysis/physics/heat.steadystate.nonlinear.h"
#include "analysis/physics/heat.transient.linear.h"
#include "analysis/physics/structuralmechanics.harmonic.real.linear.h"
#include "analysis/physics/structuralmechanics.steadystate.linear.h"
#include "analysis/physics/structuralmechanics.steadystate.nonlinear.h"

#include <fstream>

namespace espreso {

template <typename T>
FETILinearSystemSolver<T>::FETILinearSystemSolver(PhysicsConfiguration &physics, LoadStepSolverConfiguration &loadStep)
: physics(physics), loadStep(loadStep), feti(loadStep.feti)
{
    LinearSystemSolver<T>::A = &A;
    LinearSystemSolver<T>::x = &x;
    LinearSystemSolver<T>::b = &b;
    LinearSystemSolver<T>::dirichlet = &dirichlet;
}

template <typename T>
FETILinearSystemSolver<T>::~FETILinearSystemSolver()
{

}

template <typename T>
void FETILinearSystemSolver<T>::set(step::Step &step)
{
    eslog::startln("FETI: SETTING LINEAR SYSTEM", "FETI[SET]");
    feti.K.resize(A.domains.size());
    feti.x.resize(A.domains.size());
    feti.f.resize(A.domains.size());
    for (size_t di = 0; di < A.domains.size(); ++di) {
        feti.K[di].shallowCopy(A.domains[di]);
        feti.f[di].shallowCopy(b.domains[di]);
        feti.x[di].shallowCopy(x.domains[di]);
    }
    feti.decomposition = A.decomposition;
    Constrains<T>::set(step, feti, dirichlet);
    eslog::checkpointln("FETI: SET B1");
    Regularization<T>::set(step, feti);
    eslog::checkpointln("FETI: SET KERNELS");
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    feti.set(step);
    eslog::endln("FETI: LINEAR SYSTEM SET");
}

template <typename T>
void FETILinearSystemSolver<T>::update(step::Step &step)
{
    eslog::startln("FETI: UPDATING LINEAR SYSTEM", "FETI[UPDATE]");
    Constrains<T>::update(step, feti, dirichlet);
    eslog::checkpointln("FETI: UPDATE B1");
    if (A.updated) {
        Regularization<T>::update(step, feti);
        eslog::checkpointln("FETI: UPDATE KERNELS");
    }
    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: system/{K, f, R, RegMat}\n");
        for (size_t d = 0; d < feti.K.size(); ++d) {
            math::store(feti.K[d], utils::filename(utils::debugDirectory(step) + "/system", "K" + std::to_string(d)).c_str());
            math::store(feti.f[d], utils::filename(utils::debugDirectory(step) + "/system", "f" + std::to_string(d)).c_str());
            math::store(feti.R1[d], utils::filename(utils::debugDirectory(step) + "/system", "R" + std::to_string(d)).c_str());
            math::store(feti.RegMat[d], utils::filename(utils::debugDirectory(step) + "/system", "RegMat" + std::to_string(d)).c_str());
        }

        eslog::storedata(" STORE: system/{B1, B1c, B1Duplication, D2C, CMAP}\n");
        math::store(feti.c, utils::filename(utils::debugDirectory(step) + "/system", "B1c").c_str());
        for (size_t d = 0; d < feti.B1.size(); ++d) {
            math::store(feti.B1[d], utils::filename(utils::debugDirectory(step) + "/system", "B1" + std::to_string(d)).c_str());
            math::store(feti.D2C[d], utils::filename(utils::debugDirectory(step) + "/system", "D2C" + std::to_string(d)).c_str());
        }
        std::ofstream os(utils::filename(utils::debugDirectory(step) + "/system", "CMAP.txt").c_str());
        for (size_t i = 0; i < feti.lambdas.cmap.size(); ) {
            os << feti.lambdas.cmap[i] << "x:";
            for (esint d = 0; d < feti.lambdas.cmap[i + 1]; ++d) {
                os << " " << feti.lambdas.cmap[i + 2 + d];
            }
            os << "\n";
            i += feti.lambdas.cmap[i + 1] + 2;
        }
    }
    feti.updated.K = A.updated;
    feti.updated.B = step.substep == 0;
    feti.update(step);
    A.updated = b.updated = dirichlet.updated = false;
    eslog::endln("FETI: LINEAR SYSTEM UPDATED");
}

template <typename T>
bool FETILinearSystemSolver<T>::solve(step::Step &step)
{
    eslog::startln("FETI: RUN LINEAR SYSTEM", "FETI[SOLVE]");
    if (feti.solve(step)) {
        if (info::ecf->output.print_matrices) {
            eslog::storedata(" STORE: system/{x}\n");
            for (size_t d = 0; d < feti.x.size(); ++d) {
                math::store(feti.x[d], utils::filename(utils::debugDirectory(step) + "/system", "x" + std::to_string(d)).c_str());
            }
        }
        eslog::endln("FETI: LINEAR SYSTEM SOLVED");
        return true;
    }
    eslog::endln("FETI: LINEAR SYSTEM SOLVED");
    return false;
}

template struct FETILinearSystemSolver<double>;

}

