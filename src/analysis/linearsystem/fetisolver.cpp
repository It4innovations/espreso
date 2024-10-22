
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
: physics(physics), loadStep(loadStep), feti(loadStep.feti), bem(false)
{
    LinearSystemSolver<T>::A = &A;
    LinearSystemSolver<T>::x = &x.physics;
    LinearSystemSolver<T>::b = &b.physics;
    LinearSystemSolver<T>::dirichlet = &dirichlet;

    for (auto disc = info::ecf->heat_transfer.discretization.cbegin(); disc != info::ecf->heat_transfer.discretization.cend(); ++disc) {
        if (disc->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
            bem = true;
        }
    }
    for (auto disc = info::ecf->structural_mechanics.discretization.cbegin(); disc != info::ecf->structural_mechanics.discretization.cend(); ++disc) {
        if (disc->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
            bem = true;
        }
    }
}

template <typename T>
FETILinearSystemSolver<T>::~FETILinearSystemSolver()
{

}

template <typename T>
T FETILinearSystemSolver<T>::rhs_norm()
{
    Vector_FETI<Vector_Dense, T> *rhs = b.feti.copyPattern();
    for (size_t di = 0; di < b.feti.domains.size(); ++di) {
        math::copy(rhs->domains[di], b.feti.domains[di]);
        math::add(rhs->domains[di], T{-1}, feti.BtL[di]);
    }
    T norm = rhs->norm();
    delete rhs;
    return norm;
}

template <typename T>
void FETILinearSystemSolver<T>::set(step::Step &step)
{
    eslog::startln("FETI: SETTING LINEAR SYSTEM", "FETI[SET]");
    feti.K.resize(A.domains.size());
    feti.x.resize(A.domains.size());
    feti.f.resize(A.domains.size());

    b.feti.decomposition = x.feti.decomposition = A.decomposition;
    b.feti.domains.resize(A.domains.size());
    x.feti.domains.resize(A.domains.size());
    for (size_t di = 0; di < A.domains.size(); ++di) {
        b.feti.domains[di].resize(A.decomposition->dsize[di]);
        x.feti.domains[di].resize(A.decomposition->dsize[di]);
        feti.K[di].shallowCopy(A.domains[di]);
        feti.f[di].shallowCopy(b.feti.domains[di]);
        feti.x[di].shallowCopy(x.feti.domains[di]);
    }
    feti.decomposition = A.decomposition;
    regularization.set(step, feti);
    eslog::checkpointln("FETI: SET KERNELS");
    constrains.set(step, feti, dirichlet);
    eslog::checkpointln("FETI: SET B1");
    eslog::info(" = ----------------------------------------------------------------------------------------- = \n");
    feti.set(step);
    eslog::endln("FETI: LINEAR SYSTEM SET");
}

template <typename T>
void FETILinearSystemSolver<T>::update(step::Step &step)
{
    eslog::startln("FETI: UPDATING LINEAR SYSTEM", "FETI[UPDATE]");
    this->b.physics.spliteTo(&this->b.feti);

    regularization.update(step, feti);
    eslog::checkpointln("FETI: UPDATE KERNELS");

    constrains.update(step, feti, dirichlet);
    eslog::checkpointln("FETI: UPDATE B1");

    if (bem) {
        #pragma omp parallel for
        for (size_t d = 0; d < feti.K.size(); ++d) {
            Matrix_Dense<T> R, RRt, P, K, KP;
            R.resize(feti.R1[d].nrows, feti.R1[d].ncols);
            RRt.resize(R.ncols, R.ncols);
            math::copy(R, feti.R1[d]);
            math::orthonormalize(R);
            math::blas::AAt(R, RRt, true);
            for (int r = 0; r < RRt.nrows; ++r) { // make R full symmetric matrix
                for (int c = r + 1; c < RRt.ncols; ++c) {
                    RRt.vals[c * RRt.ncols + r] = RRt.vals[r * RRt.ncols + c];
                }
            }
            P.resize(feti.K[d].nrows, feti.K[d].nrows);
            KP.resize(feti.K[d].nrows, feti.K[d].nrows);
            math::eye(P, 1.);
            math::add(P, -1., RRt);
            math::copy(K, feti.K[d]);
            math::blas::multiply(T{1.}, K, P, T{0}, KP);
            math::blas::multiply(T{1.}, P, KP, T{0}, K);
            math::copy(feti.K[d], K);
        }
        eslog::checkpointln("FETI: K PROJECTED");
    }

    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: system/{K, f, R, RegMat}\n");
        if (feti.K.size() == 1) {
            math::store(feti.K[0], utils::filename(utils::debugDirectory(step) + "/system", "K").c_str());
            math::store(feti.f[0], utils::filename(utils::debugDirectory(step) + "/system", "f").c_str());
        } else {
            for (size_t d = 0; d < feti.K.size(); ++d) {
                math::store(feti.K[d], utils::filename(utils::debugDirectory(step) + "/system", "K" + std::to_string(d)).c_str());
                math::store(feti.f[d], utils::filename(utils::debugDirectory(step) + "/system", "f" + std::to_string(d)).c_str());
            }
        }
        for (size_t d = 0; d < feti.K.size(); ++d) {
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

        if (feti.B0.size()) {
            eslog::storedata(" STORE: system/{B0, D2C0}\n");
            for (size_t d = 0; d < feti.B0.size(); ++d) {
                math::store(feti.B0[d], utils::filename(utils::debugDirectory(step) + "/system", "B0" + std::to_string(d)).c_str());
                math::store(feti.D2C0[d], utils::filename(utils::debugDirectory(step) + "/system", "D2C0" + std::to_string(d)).c_str());
            }
        }
    }
    feti.updated.K = A.updated;
    feti.updated.B |= step.substep == 0;
    feti.update(step);
    A.updated = b.physics.updated = b.feti.updated = dirichlet.updated = false;
    feti.updated.K = feti.updated.B = false;
    eslog::endln("FETI: LINEAR SYSTEM UPDATED");
}

template <typename T>
bool FETILinearSystemSolver<T>::solve(step::Step &step)
{
    eslog::startln("FETI: RUN LINEAR SYSTEM", "FETI[SOLVE]");

    if (false) {
        printf("RHS\n");
        auto dd = info::mesh->nodes->domains->begin();
        for (esint n = 0; n < info::mesh->nodes->size; ++n, ++dd) {
            std::vector<Point> pp(dd->size());
            for (int d = 0; d < 3; ++d) {
                auto dmap = feti.decomposition->dmap->cbegin() + 3 * n + d;
                int ii = 0;
                for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                    if (feti.decomposition->ismy(di->domain)) {
                        pp[ii++][d] = feti.f[di->domain - feti.decomposition->dbegin].vals[di->index];
                    }
                }
            }
            for (size_t i = 1; i < pp.size(); ++i) {
                pp[0] += pp[i];
            }
            printf("%2d [%+.14e %+.14e %+.14e]\n", n, pp[0].x, pp[0].y, pp[0].z);
        }
    }

    bool result = feti.solve(step);
    constrains.eq.enforce(step, feti, dirichlet);

    if (true) {
        double x1 = 1e12, x2 = 1 / x1;
        for (size_t d = 0; d < feti.x.size(); ++d) {
            for (int i = 0; i < feti.x[d].size; ++i) {
                if (std::fabs(feti.x[d].vals[i]) < x2) {
                    feti.x[d].vals[i] = 0;
                } else {
                    feti.x[d].vals[i] = std::ceil(x1 * feti.x[d].vals[i]) * x2;
                }
            }
        }
    }

    this->x.feti.averageTo(&this->x.physics);

    if (info::ecf->output.print_matrices) {
        eslog::storedata(" STORE: system/{x}\n");
        if (feti.x.size() == 1) {
            math::store(feti.x[0], utils::filename(utils::debugDirectory(step) + "/system", "x").c_str());
        } else {
            for (size_t d = 0; d < feti.x.size(); ++d) {
                math::store(feti.x[d], utils::filename(utils::debugDirectory(step) + "/system", "x" + std::to_string(d)).c_str());
            }
        }
    }

    if (false) {
        printf("SOLUTION\n");
        auto dd = info::mesh->nodes->domains->begin();
        for (esint n = 0; n < info::mesh->nodes->size; ++n, ++dd) {
            std::vector<Point> pp(dd->size());
            for (int d = 0; d < 3; ++d) {
                auto dmap = feti.decomposition->dmap->cbegin() + 3 * n + d;
                int ii = 0;
                for (auto di = dmap->begin(); di != dmap->end(); ++di) {
                    if (feti.decomposition->ismy(di->domain)) {
                        pp[ii++][d] = feti.x[di->domain - feti.decomposition->dbegin].vals[di->index];
                    }
                }
            }
            for (size_t i = 0; i < pp.size(); ++i) {
                printf("%2d [%+.14e %+.14e %+.14e]\n", n, pp[i].x, pp[i].y, pp[i].z);
            }
        }
    }

    eslog::endln("FETI: LINEAR SYSTEM SOLVED");
    return result;
}

template struct FETILinearSystemSolver<double>;

}

