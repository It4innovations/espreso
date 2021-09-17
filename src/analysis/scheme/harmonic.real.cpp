
#include "analysis/scheme/harmonic.real.h"
#include "analysis/linearsystem/linearsystem.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/physicssolver/harmonic.h"

#include <cmath>

using namespace espreso;


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

	system->solver.A->add(-frequency.angular, C, 0, dofs, dofs, 2 * dofs);
	system->solver.A->add( frequency.angular, C, dofs, 0, dofs, 2 * dofs);

	system->solver.b->touched = true;
	system->solver.b->fill(0);
	system->solver.b->add(1., re.f, 0   , dofs, 2 * dofs);
	system->solver.b->add(1., im.f, dofs, dofs, 2 * dofs);

	system->solver.dirichlet->touched = true;
	system->solver.dirichlet->fill(0);
	system->solver.dirichlet->add(1, re.dirichlet, 0   , dofs, 2 * dofs);
	system->solver.dirichlet->add(1, im.dirichlet, dofs, dofs, 2 * dofs);

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{K, M, C, f.re, f.im, dirichlet.re, dirichlet.im}\n");
		K->store(utils::filename(utils::debugDirectory() + "/scheme", "K").c_str());
		M->store(utils::filename(utils::debugDirectory() + "/scheme", "M").c_str());
		C->store(utils::filename(utils::debugDirectory() + "/scheme", "C").c_str());
		re.f->store(utils::filename(utils::debugDirectory() + "/scheme", "f.re").c_str());
		im.f->store(utils::filename(utils::debugDirectory() + "/scheme", "f.im").c_str());
		re.dirichlet->store(utils::filename(utils::debugDirectory() + "/scheme", "dirichlet.re").c_str());
		im.dirichlet->store(utils::filename(utils::debugDirectory() + "/scheme", "dirichlet.im").c_str());
	}
}

void AX_HarmonicReal::extractSolution(step::Frequency &frequency, AX_LinearSystem<double> *system)
{
	re.x->fillData(system->solver.x, 0   , dofs, 2 * dofs);
	im.x->fillData(system->solver.x, dofs, dofs, 2 * dofs);

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{x.re, x.im}\n");
		re.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.re").c_str());
		im.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.im").c_str());
	}
}

