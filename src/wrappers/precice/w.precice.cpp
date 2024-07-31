
#include "w.precice.h"

#include "esinfo/mpiinfo.h"

using namespace espreso;

#ifdef HAVE_PRECICE
#include <precice/precice.hpp>
#endif

Precice::Precice()
{
#ifdef HAVE_PRECICE
    precice::Participant precice("SolidSolver", "precice-config.xml", info::mpi::rank, info::mpi::size); // constructor
#endif
}

Precice::~Precice()
{
#ifdef HAVE_PRECICE

#endif
}
