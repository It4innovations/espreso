
#include "harmonic.complex.h"
#include "analysis/linearsystem/linearsystem.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"

#include <cmath>

using namespace espreso;

void AX_HarmonicComplex::composeSystem(step::Frequency &frequency, AX_LinearSystem<double, std::complex<double> > *system)
{
	#if 0
	// A = [
	// K - omega^2 * M                -iC
	//              iC    K - omega^2 * M
	// ]
	system->solver.A->touched = true;
	system->solver.A->fill(0);
	system->solver.A->sum(1., K, -frequency.angular * frequency.angular, M, 0   , 0   , dofs, 2 * dofs);
	system->solver.A->sum(1., K, -frequency.angular * frequency.angular, M, dofs, dofs, dofs, 2 * dofs);

	system->solver.A->add(-frequency.angular, C, 0, dofs, dofs, 2 * dofs);
	system->solver.A->add(frequency.angular, C, dofs, 0, dofs, 2 * dofs);

	system->solver.b->touched = true;
	system->solver.b->fill(0);
	system->solver.b->add(1., re.f, 0   , dofs, 2 * dofs);
	system->solver.b->add(1., im.f, dofs, dofs, 2 * dofs);
	#endif

	// system->solver.A->touched = true;
	// system->solver.A->fill(std::complex<double>(0.0));
	// system->solver.A->sum(std::complex<double>(1.0), ???)
}

void AX_HarmonicComplex::extractSolution(step::Frequency &frequency, AX_LinearSystem<double, std::complex<double> > *system)
{

}

