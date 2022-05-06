
#include "scheme.h"
#include "harmonic.h"
#include "analysis/linearsystem/linearsystem.h"
#include "esinfo/eslog.h"
#include "esinfo/stepinfo.h"
#include "config/ecf/physics/physicssolver/harmonic.h"

#include <cmath>

using namespace espreso;

Harmonic::Harmonic(HarmonicSolverConfiguration &configuration, int dofs)
: configuration(configuration), dofs(dofs), K{}, M{}, C{}, re{}, im{}
{

}



Harmonic::~Harmonic()
{
	clear(K, M, C, re.f, im.f, re.x, im.x, re.dirichlet, im.dirichlet);
}

void Harmonic::initFrequency(step::Frequency &frequency)
{
	frequency.shift = (configuration.max_frequency - configuration.min_frequency) / configuration.num_samples;
	frequency.start = configuration.min_frequency;
	frequency.current = frequency.start;
	frequency.angular = 2 * M_PI * frequency.current;
	frequency.final = configuration.max_frequency;
}

void Harmonic::nextFrequency(step::Frequency &frequency)
{
	frequency.current += frequency.shift;
	if (frequency.current + frequency.precision >= frequency.final) {
		frequency.current = frequency.final;
	}
	frequency.angular = 2 * M_PI * frequency.current;
}
