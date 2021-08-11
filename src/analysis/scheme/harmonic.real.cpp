
#include "analysis/scheme/harmonic.real.h"
#include "analysis/linearsystem/linearsystem.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/physicssolver/harmonic.h"

#include <cmath>

using namespace espreso;

AX_HarmonicReal::AX_HarmonicReal(HarmonicSolverConfiguration &configuration, int dofs)
: configuration(configuration), dofs(dofs), K{}, M{}, C{}, re{}, im{}
{

}

AX_HarmonicReal::~AX_HarmonicReal()
{
	if (K) {
		delete K;
	}
	if (M) {
		delete M;
	}
	if (C) {
		delete C;
	}
	if (re.f) {
		delete re.f;
	}
	if (re.x) {
		delete re.x;
	}
	if (im.f) {
		delete im.f;
	}
	if (im.x) {
		delete im.x;
	}
}

void AX_HarmonicReal::initFrequency(step::Frequency &frequency)
{
	frequency.shift = (configuration.max_frequency - configuration.min_frequency) / configuration.num_samples;
	frequency.start = configuration.min_frequency;
	frequency.current = frequency.start;
	frequency.angular = 2 * M_PI * frequency.current;
	frequency.final = configuration.max_frequency;
}

void AX_HarmonicReal::nextFrequency(step::Frequency &frequency)
{
	frequency.current += frequency.shift;
	if (frequency.current + frequency.precision >= frequency.final) {
		frequency.current = frequency.final;
	}
	frequency.angular = 2 * M_PI * frequency.current;
}

void AX_HarmonicReal::init(AX_LinearSystem<double> *system)
{
	system->setMapping(K = system->assembler.A->copyPattern());
	system->setMapping(M = system->assembler.A->copyPattern());
//	system->setMapping(C = system->assembler.A->copyPattern());
	system->setMapping(re.f = system->assembler.b->copyPattern());
	system->setMapping(re.x = system->assembler.b->copyPattern());
	system->setMapping(im.f = system->assembler.b->copyPattern());
	system->setMapping(im.x = system->assembler.b->copyPattern());
}

void AX_HarmonicReal::composeSystem(step::Frequency &frequency, AX_LinearSystem<double> *system)
{
	// A = [
	// K - omega^2 * M                -iC
	//              iC    K - omega^2 * M
	// ]
	system->solver.A->touched = true;
	system->solver.A->fill(0);
	system->solver.A->sum(1., K, -frequency.angular * frequency.angular, M, 0   , 0   , dofs, 2 * dofs);
	system->solver.A->sum(1., K, -frequency.angular * frequency.angular, M, dofs, dofs, dofs, 2 * dofs);

	system->solver.b->touched = true;
	system->solver.b->fill(0);
	system->solver.b->add(1., re.f, 0   , dofs, 2 * dofs);
	system->solver.b->add(1., im.f, dofs, dofs, 2 * dofs);
}

void AX_HarmonicReal::composeDirichlet(AX_LinearSystem<double> *system)
{
	math::fill(*system->solver.dirichlet, 0.);
	math::add(*system->solver.dirichlet, 1., *system->assembler.dirichlet, 0   , dofs, 2 * dofs);
	math::add(*system->solver.dirichlet, 1., *system->assembler.dirichlet, dofs, dofs, 2 * dofs);
}

void AX_HarmonicReal::extractSolution(AX_LinearSystem<double> *system)
{
	re.x->fillData(system->solver.x, 0   , dofs, 2 * dofs);
	im.x->fillData(system->solver.x, dofs, dofs, 2 * dofs);
}

void AX_HarmonicReal::storeScheme(step::Frequency &frequency)
{
	if (info::ecf->output.print_matrices) {
		K->store(utils::filename(utils::debugDirectory() + "/scheme", "K").c_str());
		M->store(utils::filename(utils::debugDirectory() + "/scheme", "M").c_str());
		re.f->store(utils::filename(utils::debugDirectory() + "/scheme", "f.re").c_str());
		im.f->store(utils::filename(utils::debugDirectory() + "/scheme", "f.im").c_str());
	}
}

void AX_HarmonicReal::storeSolution(step::Frequency &frequency)
{
	if (info::ecf->output.print_matrices) {
		re.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.re").c_str());
		im.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.im").c_str());
	}
}

