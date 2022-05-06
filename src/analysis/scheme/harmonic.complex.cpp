
#include "harmonic.complex.h"
#include "analysis/linearsystem/linearsystem.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/stepinfo.h"

#include <cmath>

using namespace espreso;

void HarmonicComplex::composeSystem(step::Frequency &frequency, LinearSystem<double, std::complex<double> > *system)
{
	// A = K - omega^2 * M + iC
	system->solver.A->touched = true;
	system->solver.A->set(std::complex<double>(0, 0));
	system->solver.A->copyReal(K);
	system->solver.A->addReal(-frequency.angular * frequency.angular, M);
	system->solver.A->addImag(frequency.angular, C);

	system->solver.b->touched = true;
	system->solver.b->copyReal(re.f);
	system->solver.b->copyImag(im.f);

	system->solver.dirichlet->touched = true;
	system->solver.dirichlet->set(std::complex<double>(0, 0));
	system->solver.dirichlet->copyReal(re.dirichlet);
	system->solver.dirichlet->copyImag(im.dirichlet);

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

void HarmonicComplex::extractSolution(step::Frequency &frequency, LinearSystem<double, std::complex<double> > *system)
{
	system->solver.x->copyRealTo(re.x);
	system->solver.x->copyImagTo(im.x);

	if (info::ecf->output.print_matrices) {
		eslog::storedata(" STORE: scheme/{x.re, x.im}\n");
		re.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.re").c_str());
		im.x->store(utils::filename(utils::debugDirectory() + "/scheme", "x.im").c_str());
	}
}


