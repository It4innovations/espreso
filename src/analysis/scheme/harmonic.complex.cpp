
#include "harmonic.complex.h"
#include "analysis/linearsystem/linearsystem.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"

#include <cmath>

using namespace espreso;

void AX_HarmonicComplex::composeSystem(step::Frequency &frequency, AX_LinearSystem<double, std::complex<double> > *system)
{
	// A = K - omega^2 * M + iC
	system->solver.A->touched = true;
	system->solver.A->fill(0);

	system->solver.A->sum(1., K, -frequency.angular * frequency.angular, M);
	system->solver.A->add_imag(frequency.angular, C);

	system->solver.b->touched = true;
	system->solver.b->fill(0);
	system->solver.b->add(1., re.f);
	system->solver.b->add_imag(1., im.f);
}

void AX_HarmonicComplex::extractSolution(step::Frequency &frequency, AX_LinearSystem<double, std::complex<double> > *system)
{
	math::fill(*system->solver.dirichlet, std::complex<double>(0.0));
	math::add(*system->solver.dirichlet, 1.0, *system->assembler.dirichlet);
}

}

