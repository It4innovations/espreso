
#include "w.precice.h"

#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/eslog.h"
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
    size_t size;
    esint *ids;
    std::vector<double> data;
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

        _data->size = info::mesh->surface->nIDs->datatarray().size();
        _data->ids = info::mesh->surface->nIDs->datatarray().data();

        size_t size = info::mesh->surface->nIDs->datatarray().size();
        if (info::mesh->dimension == 2) {
            std::vector<double> coords; coords.reserve(info::mesh->dimension * size);
            for (size_t n = 0; n < size; ++n) {
                coords.push_back(info::mesh->surface->coordinates->datatarray()[n].x);
                coords.push_back(info::mesh->surface->coordinates->datatarray()[n].y);
            }
            _data->precice.setMeshVertices("SolidMesh", coords, precice::span(_data->ids, _data->ids + _data->size));
        } else {
            double *coords = &info::mesh->surface->coordinates->datatarray().data()->x;
            _data->precice.setMeshVertices("SolidMesh", precice::span(coords, coords + info::mesh->dimension * _data->size), precice::span(_data->ids, _data->ids + _data->size));
        }
        _data->precice.initialize();
        _data->data.resize(info::mesh->dimension * _data->size);
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
        _data->precice.readData("SolidMesh", name, precice::span(_data->ids, _data->ids + _data->size), dt, _data->data);
        for (size_t n = 0; n < _data->size; ++n) {
            for (int d = 0; d < info::mesh->dimension; ++d) {
                data[_data->ids[n] * info::mesh->dimension + d] = _data->data[n * info::mesh->dimension + d];
            }
        }
    }
#endif
}

void Precice::write(const char *name, double *data)
{
#ifdef HAVE_PRECICE
    if (_data) {
        for (size_t n = 0; n < _data->size; ++n) {
            for (int d = 0; d < info::mesh->dimension; ++d) {
                _data->data[n * info::mesh->dimension + d] = data[_data->ids[n] * info::mesh->dimension + d];
            }
        }
        _data->precice.writeData("SolidMesh", name, precice::span(_data->ids, _data->ids + _data->size), _data->data);
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

    size_t size = info::mesh->surface->nIDs->datatarray().size();
    esint *ids = info::mesh->surface->nIDs->datatarray().data();
    if (info::mesh->dimension == 2) {
        std::vector<double> coords; coords.reserve(info::mesh->dimension * size);
        for (size_t n = 0; n < size; ++n) {
            coords.push_back(info::mesh->surface->coordinates->datatarray()[n].x);
            coords.push_back(info::mesh->surface->coordinates->datatarray()[n].y);
        }
        participant.setMeshVertices("FluidMesh", coords, precice::span(ids, ids + size));
    } else {
        double *coords = &info::mesh->surface->coordinates->datatarray().data()->x;
        participant.setMeshVertices("FluidMesh", precice::span(coords, coords + info::mesh->dimension * size), precice::span(ids, ids + size));
    }

    std::vector<double> readData(size * info::mesh->dimension);
    std::vector<double> writeData(size * info::mesh->dimension);

    for (size_t i = info::mesh->dimension - 1; i < size * info::mesh->dimension; i += info::mesh->dimension) {
        writeData[i] = 1;
    }
    if (participant.requiresInitialData()) {
        participant.writeData(meshName, dataWriteName, precice::span(ids, ids + size), writeData);
    }
    participant.initialize();

    while (participant.isCouplingOngoing()) {
        if (participant.requiresWritingCheckpoint()) {
            std::cout << "DUMMY: Writing iteration checkpoint\n";
        }

        double dt = participant.getMaxTimeStepSize();
        participant.readData(meshName, dataReadName, precice::span(ids, ids + size), dt, readData);

        for (size_t i = info::mesh->dimension - 1; i < size * info::mesh->dimension; i += info::mesh->dimension) {
            writeData[i] = readData[i];
        }

        participant.writeData(meshName, dataWriteName, precice::span(ids, ids + size), writeData);

//        printf("DISPLACEMENT\n");
//        for (int d = 0; d < info::mesh->dimension; ++d) {
//            for (size_t i = d; i < size * info::mesh->dimension; i += info::mesh->dimension) {
//                printf(" %+e", readData[i]);
//            }
//            printf("\n");
//        }
//
//        printf("FORCES\n");
//        for (int d = 0; d < info::mesh->dimension; ++d) {
//            for (size_t i = d; i < size * info::mesh->dimension; i += info::mesh->dimension) {
//                printf(" %+e", writeData[i]);
//            }
//            printf("\n");
//        }

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
