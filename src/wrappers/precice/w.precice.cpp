
#include "w.precice.h"

#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/surfacestore.h"

#include <iostream>

#ifdef HAVE_PRECICE
#include <precice/precice.hpp>

namespace espreso {
struct PreciceData {
    PreciceData(const char *name)
    : precice(name, "precice-config.xml", info::mpi::rank, info::mpi::size)
    {

    }

    precice::Participant precice;
};
}

#endif

using namespace espreso;

Precice::Precice(const char *name, bool active)
: _data(nullptr)
{
#ifdef HAVE_PRECICE
    if (active) {
        _data = new PreciceData(name);

        size_t size = info::mesh->surface->size;
        double *coords = &info::mesh->surface->coordinates->datatarray().data()->x;
        esint *ids = info::mesh->surface->nIDs->datatarray().data();
        _data->precice.setMeshVertices("SolidMesh", precice::span(coords, coords + 3 * size), precice::span(ids, ids + size));
    }
#endif
}

Precice::~Precice()
{
#ifdef HAVE_PRECICE
    if (_data) {
        delete _data;
    }
#endif
}

void Precice::read(const char *name, double *data, double dt)
{
#ifdef HAVE_PRECICE
    if (_data) {
        size_t size = info::mesh->surface->size;
        esint *ids = info::mesh->surface->nIDs->datatarray().data();
        _data->precice.readData("SolidMesh", name, precice::span(ids, ids + size), dt, precice::span(data, data + 3 * size));
    }
#endif
}

void Precice::write(const char *name, double *data)
{
#ifdef HAVE_PRECICE
    if (_data) {
        size_t size = info::mesh->surface->size;
        esint *ids = info::mesh->surface->nIDs->datatarray().data();
        _data->precice.writeData("SolidMesh", name, precice::span(ids, ids + size), precice::span(data, data + 3 * size));
    }
#endif
}

void Precice::advance(double dt)
{
#ifdef HAVE_PRECICE
    if (_data) {
        _data->precice.advance(dt);
    }
#endif
}

void Precice::dummy()
{
#ifdef HAVE_PRECICE
    using namespace precice;

    std::string configFileName("precice-config.xml");
    std::string solverName("FluidSolver");
    std::string meshName("FluidMesh");
    std::string dataWriteName("Forces");
    std::string dataReadName("Displacement");

    std::cout << "DUMMY: Running solver dummy with preCICE config file \"" << configFileName << "\" and participant name \"" << solverName << "\".\n";

    Participant participant(solverName, configFileName, info::mpi::rank, info::mpi::size);

    size_t size = info::mesh->surface->size;
    double *coords = &info::mesh->surface->coordinates->datatarray().data()->x;
    esint *ids = info::mesh->surface->nIDs->datatarray().data();
    participant.setMeshVertices("SolidMesh", precice::span(coords, coords + 3 * size), precice::span(ids, ids + size));

    std::vector<double> readData(size * 3);
    std::vector<double> writeData(size * 3);

    participant.initialize();

    while (participant.isCouplingOngoing()) {
      if (participant.requiresWritingCheckpoint()) {
        std::cout << "DUMMY: Writing iteration checkpoint\n";
      }

      double dt = participant.getMaxTimeStepSize();
      participant.readData(meshName, dataReadName, precice::span(ids, ids + size), dt, readData);

      for (size_t i = 2; i < size * 3; i += 3) {
        writeData.at(i) = readData.at(i) * 1.5;
      }

      participant.writeData(meshName, dataWriteName, precice::span(ids, ids + size), writeData);

      participant.advance(dt);

      if (participant.requiresReadingCheckpoint()) {
        std::cout << "DUMMY: Reading iteration checkpoint\n";
      } else {
        std::cout << "DUMMY: Advancing in time\n";
      }
    }

    participant.finalize();
    std::cout << "DUMMY: Closing C++ solver dummy...\n";
#endif
}
