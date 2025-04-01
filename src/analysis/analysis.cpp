
#include "analysis.h"
#include "physics/acoustic.real.linear.h"
#include "physics/acoustic.complex.linear.h"
#include "physics/heat.steadystate.linear.h"
#include "physics/heat.steadystate.nonlinear.h"
#include "physics/heat.transient.linear.h"
#include "physics/structuralmechanics.steadystate.linear.h"
#include "physics/structuralmechanics.steadystate.nonlinear.h"
#include "physics/structuralmechanics.transient.linear.h"
#include "physics/structuralmechanics.transient.nonlinear.h"
#include "physics/structuralmechanics.harmonic.real.linear.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "esinfo/stepinfo.h"
#include "output/output.h"

using namespace espreso;

void Analysis::run()
{
    eslog::startln("ESPRESO: SIMULATION STARTED", "SIMULATION");

//    for (auto range = info::ecf->ranges.begin(); range != info::ecf->ranges.end(); ++range) {
//
//    }

    step::Step step;
    Physics *current = nullptr, *prev = nullptr;

    auto solve = [&] () {
        if (current == nullptr) {
            eslog::globalerror("not implemented physics\n");
        }

        if (!current->analyze(step)) {
            eslog::globalerror("physical analysis failed\n");
        }
        eslog::checkpointln("SIMULATION: PHYSICS ANALYSED");
        if (!current->run(step, prev)) {
            #warning "temporary commented out for benchmarks, uncomment for production"
            // eslog::globalerror("physical solver failed\n");
        }
        eslog::checkpointln("SIMULATION: PHYSICS SOLVED");

        if (prev) delete prev;
        prev = current;
        current = nullptr;
    };

    switch (info::ecf->physics) {
//    case PhysicsConfiguration::TYPE::ACOUSTICS:
//        switch (info::ecf->acoustics.load_steps_settings.at(s).system) {
//        case AcousticLoadStepConfiguration::SYSTEM::REAL: physics = new AcousticRealLinear(info::ecf->acoustics, info::ecf->acoustics.load_steps_settings.at(s)); break;
//        case AcousticLoadStepConfiguration::SYSTEM::COMPLEX: physics = new AcousticComplexLinear(info::ecf->acoustics, info::ecf->acoustics.load_steps_settings.at(s)); break;
//        }
//        break;
    case PhysicsConfiguration::TYPE::HEAT_TRANSFER:
        step.loadsteps = info::ecf->heat_transfer.load_steps;
        for (int s = 1; s <= step.loadsteps; ++s, ++step.loadstep) {
            switch (info::ecf->heat_transfer.load_steps_settings.at(s).type) {
            case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
                switch (info::ecf->heat_transfer.load_steps_settings.at(s).mode) {
                case LoadStepSolverConfiguration::MODE::LINEAR: current = new HeatSteadyStateLinear(info::ecf->heat_transfer, info::ecf->heat_transfer.load_steps_settings.at(s)); break;
                case LoadStepSolverConfiguration::MODE::NONLINEAR: current = new HeatSteadyStateNonLinear(info::ecf->heat_transfer, info::ecf->heat_transfer.load_steps_settings.at(s)); break;
                } break;
            case LoadStepSolverConfiguration::TYPE::TRANSIENT:
                switch (info::ecf->heat_transfer.load_steps_settings.at(s).mode) {
                case LoadStepSolverConfiguration::MODE::LINEAR: current = new HeatTransientLinear(info::ecf->heat_transfer, info::ecf->heat_transfer.load_steps_settings.at(s)); break;
                case LoadStepSolverConfiguration::MODE::NONLINEAR: eslog::globalerror("implement HeatTransientNonLinear solver.\n"); break; // physics = new HeatSteadyStateNonLinear(info::ecf->heat_transfer, info::ecf->heat_transfer.load_steps_settings.at(s)); break;
                } break;
            case LoadStepSolverConfiguration::TYPE::HARMONIC:
                eslog::globalerror("invalid combination: HARMONIC -- HEAT_TRANSFER.\n");
                break;
            }

            solve();
        }
        break;
    case PhysicsConfiguration::TYPE::STRUCTURAL_MECHANICS:
        step.loadsteps = info::ecf->structural_mechanics.load_steps;
        for (int s = 1; s <= step.loadsteps; ++s, ++step.loadstep) {
            switch (info::ecf->structural_mechanics.load_steps_settings.at(s).type) {
            case LoadStepSolverConfiguration::TYPE::HARMONIC:
                switch (info::ecf->structural_mechanics.load_steps_settings.at(s).mode) {
                case LoadStepSolverConfiguration::MODE::LINEAR: current = new StructuralMechanicsHarmonicRealLinear(info::ecf->structural_mechanics, info::ecf->structural_mechanics.load_steps_settings.at(s)); break;
                case LoadStepSolverConfiguration::MODE::NONLINEAR:  break;
                } break;
            case LoadStepSolverConfiguration::TYPE::STEADY_STATE:
                switch (info::ecf->structural_mechanics.load_steps_settings.at(s).mode) {
                case LoadStepSolverConfiguration::MODE::LINEAR: current = new StructuralMechanicsSteadyStateLinear(info::ecf->structural_mechanics, info::ecf->structural_mechanics.load_steps_settings.at(s)); break;
                case LoadStepSolverConfiguration::MODE::NONLINEAR: current = new StructuralMechanicsSteadyStateNonLinear(info::ecf->structural_mechanics, info::ecf->structural_mechanics.load_steps_settings.at(s)); break;
                } break;
            case LoadStepSolverConfiguration::TYPE::TRANSIENT:
                switch (info::ecf->structural_mechanics.load_steps_settings.at(s).mode) {
                case LoadStepSolverConfiguration::MODE::LINEAR:    current = new StructuralMechanicsTransientLinear(info::ecf->structural_mechanics, info::ecf->structural_mechanics.load_steps_settings.at(s)); break;
                case LoadStepSolverConfiguration::MODE::NONLINEAR: current = new StructuralMechanicsTransientNonLinear(info::ecf->structural_mechanics, info::ecf->structural_mechanics.load_steps_settings.at(s)); break;
                } break;
            }

            solve();
        }
        break;
    }

    if (prev) delete prev;
    eslog::endln("SIMULATION: SIMULATION FINISHED");
}
