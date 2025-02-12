
#include "topologyoptimization.h"
#include "basis/utilities/parser.h"
#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "config/ecf/physics/physicssolver/loadstep.h"
#include "esinfo/stepinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/mesh.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "math/matrix.h"
#include "math/vector.dense.h"
#include "physics/system/linearsystem.h"
#include "physics/system/builder/builder.h"
#include "physics/substepsolver/substepsolver.h"

#include <cmath>

using namespace espreso;

TopologyOptimization::TopologyOptimization(LinearSystem *system, SubStepSolver *subStepSolver, TopologyOptimizationConfiguration &configuration)
: LoadStepSolver(system, subStepSolver, 1), _configuration(configuration),
  xPhys(NULL), x(NULL), DC(NULL), C(NULL), DV(NULL)
{

}

TopologyOptimization::~TopologyOptimization()
{
    if (xPhys) { delete xPhys; }
    if (x) { delete x; }
    if (DC) { delete DC; }
    if (C) { delete C; }
    if (DV) { delete DV; }
}

void TopologyOptimization::init(LoadStepSolver *previous)
{
    for (size_t i = 0; i < info::mesh->elements->data.size(); i++) {
        if (StringCompare::caseInsensitiveEq(info::mesh->elements->data[i]->name, "DESIGN_VARIABLE")) {
            xPhys = new VectorDense(info::mesh->elements->data[i]->data.size(), info::mesh->elements->data[i]->data.data());
        }
        if (StringCompare::caseInsensitiveEq(info::mesh->elements->data[i]->name, "COMPLIANCE_DERIVATION")) {
            DC = new VectorDense(info::mesh->elements->data[i]->data.size(), info::mesh->elements->data[i]->data.data());
        }
        if (StringCompare::caseInsensitiveEq(info::mesh->elements->data[i]->name, "COMPLIANCE")) {
            C = new VectorDense(info::mesh->elements->data[i]->data.size(), info::mesh->elements->data[i]->data.data());
        }
    }
    if (xPhys == NULL) {
        eslog::internalFailure("topology optimization need DESIGN_VARIABLE.\n");
    }
    x = xPhys->shallowCopyStructure();
    x->fillData(xPhys);
    DV = xPhys->shallowCopyStructure();
    DV->fill(1);
}

void TopologyOptimization::updateStructuralMatrices()
{
    _system->builder->matrices &= Builder::Request::K | Builder::Request::RBCf;
    _system->assemble();
}

void TopologyOptimization::runNextSubstep()
{
    _system->builder->internalForceReduction = 1;
    _system->builder->timeIntegrationConstantK = 1;
    _system->builder->timeIntegrationConstantC = 0;
    _system->builder->timeIntegrationConstantM = 0;

    auto setFixedRegions = [&] () {
        for (auto preset = _configuration.constraint.preset_regions.begin(); preset != _configuration.constraint.preset_regions.end(); ++preset) {
            ElementsRegionStore *region = info::mesh->eregion(preset->first);
            double value = preset->second == TopologyOptimizationConstraint::PresetValues::SOLID ? 1 : 0;
            #pragma omp parallel for
            for (int t = 0; t < info::env::OMP_NUM_THREADS; t++) {
                for (auto e = region->elements->datatarray().cbegin(t); e != region->elements->datatarray().cend(t); ++e) {
                    xPhys->vals[*e] = value;
                }
            }
        }
    };

    setFixedRegions();

    int iteration = 1;
    double change = .5;
    do {
        step::time.current += 0.01;
        step::time.shift = 0.1;
        _system->nextSubstep();

        eslog::solver("\n = ================================= TOPOLOGY OPTIMIZATION ================================= =\n");
        eslog::solver(" =  LOAD STEP %2d, SUBSTEP %4d,                             ITERATION %4d, CHANGE %8.6f  =\n", step::step.loadstep + 1, step::step.substep + 1, iteration, change);
        eslog::solver(" = ----------------------------------------------------------------------------------------- =\n");

        _subStepSolver->solve(*this);
        _system->processSolution();

        double c = 0;
        for (esint i = 0; i < C->size; i++) {
            c += C->vals[i];
        }

        switch (_configuration.solver_settings.type) {
        case TopologyOptimizationSolverSettings::Type::OC: {
            double l1 = _configuration.solver_settings.oc.lower_bound;
            double l2 = _configuration.solver_settings.oc.upper_bound;
            double move = _configuration.solver_settings.oc.move;
            while ((l2 - l1) / (l2 + l1) > _configuration.solver_settings.oc.precision) {
                double lmid = .5 * (l2 + l1);
                double rsum = 0, sum;
                for (esint i = 0; i < xPhys->size; i++) {
                    xPhys->vals[i] = std::max(0.0, std::max(x->vals[i] - move, std::min(1.0, std::min(x->vals[i] + move, x->vals[i] * sqrt(-DC->vals[i] / DV->vals[i] / lmid)))));
                }
                setFixedRegions();
                for (esint i = 0; i < xPhys->size; i++) {
                    rsum += xPhys->vals[i];
                }
                Communication::allReduce(&rsum, &sum, 1, MPI_DOUBLE, MPI_SUM);
                if (sum > _configuration.constraint.value * info::mesh->elements->distribution.process.totalSize) {
                    l1 = lmid;
                } else {
                    l2 = lmid;
                }
            }
        } break;
        case TopologyOptimizationSolverSettings::Type::MMA: {
            eslog::internalFailure("not implemented optimization type.\n");
        } break;
        }


        double rchange = 0;
        for (esint i = 0; i < xPhys->size; i++) {
            rchange = std::max(rchange, std::fabs(xPhys->vals[i] - x->vals[i]));
        }
        Communication::allReduce(&rchange, &change, 1, MPI_DOUBLE, MPI_MAX);
        x->fillData(xPhys);

        eslog::checkpointln("PHYSICS SOLVER: DESIGN VARIABLE UPDATED");

        eslog::solver(" = ========================================================================================= =\n");
        eslog::solver(" = ================================================================= run time %12.3f s =\n\n", eslog::duration());
    } while (iteration++ < _configuration.solver_settings.max_iterations && change > _configuration.solver_settings.precision);

    step::time.current = step::time.final;
}

