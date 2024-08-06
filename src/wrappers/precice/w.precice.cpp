
#include "w.precice.h"

#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/surfacestore.h"

#include <iostream>

#ifdef HAVE_PRECICE
#include <precice/precice.hpp>

namespace espreso {
struct PreciceData {
    PreciceData()
    : precice(info::ecf->coupling.solver, info::ecf->coupling.configuration, info::mpi::rank, info::mpi::size),
      size(0)
    {

    }

    precice::Participant precice;
    size_t size;
    std::vector<esint> ids;
    std::vector<double> data;
};
}

#endif

using namespace espreso;

Precice::Precice()
: _data(nullptr)
{
#ifdef HAVE_PRECICE
    if (info::ecf->coupling.active) {
        _data = new PreciceData();

        _data->size = info::mesh->surface->nIDs->datatarray().size();
        _data->ids.insert(_data->ids.end(), info::mesh->surface->nIDs->datatarray().begin(), info::mesh->surface->nIDs->datatarray().end());

        if (info::mesh->dimension == 2) {
            std::vector<double> coords; coords.reserve(info::mesh->dimension * _data->size);
            for (size_t n = 0; n < _data->size; ++n) {
                coords.push_back(info::mesh->surface->coordinates->datatarray()[n].x);
                coords.push_back(info::mesh->surface->coordinates->datatarray()[n].y);
            }
            _data->precice.setMeshVertices(info::ecf->coupling.mesh, coords, _data->ids);
        } else {
            double *coords = &info::mesh->surface->coordinates->datatarray().data()->x;
            _data->precice.setMeshVertices(info::ecf->coupling.mesh, precice::span(coords, coords + info::mesh->dimension * _data->size), _data->ids);
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

double Precice::timeStep(double dt)
{
#ifdef HAVE_PRECICE
    if (_data) {
        return std::min(dt, _data->precice.getMaxTimeStepSize());
    }
#endif
    return dt;
}

bool Precice::requiresWritingCheckpoint()
{
#ifdef HAVE_PRECICE
    if (_data) {
        return _data->precice.requiresWritingCheckpoint();
    }
#endif
    return false;
}

bool Precice::requiresReadingCheckpoint()
{
#ifdef HAVE_PRECICE
    if (_data) {
        return _data->precice.requiresReadingCheckpoint();
    }
#endif
    return false;
}

void Precice::read(double *data, double dt)
{
#ifdef HAVE_PRECICE
    if (_data) {
        _data->precice.readData(info::ecf->coupling.mesh, info::ecf->coupling.data_in, _data->ids, dt, _data->data);
        for (size_t n = 0; n < _data->size; ++n) {
            for (int d = 0; d < info::mesh->dimension; ++d) {
                data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d] = _data->data[n * info::mesh->dimension + d];
            }
        }
    }
#endif
}

void Precice::write(double *data)
{
#ifdef HAVE_PRECICE
    if (_data) {
        for (size_t n = 0; n < _data->size; ++n) {
            for (int d = 0; d < info::mesh->dimension; ++d) {
                _data->data[n * info::mesh->dimension + d] = data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d];
            }
        }
        _data->precice.writeData(info::ecf->coupling.mesh, info::ecf->coupling.data_out, _data->ids, _data->data);
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

    std::cout << "DUMMY: Running solver dummy with preCICE config file \"" << info::ecf->coupling.dummy.configuration << "\" and participant name \"" << info::ecf->coupling.dummy.solver << "\".\n";

    Participant participant(info::ecf->coupling.dummy.solver, info::ecf->coupling.dummy.configuration, info::mpi::rank, info::mpi::size);

    size_t size = info::mesh->surface->nIDs->datatarray().size();
    esint *ids = info::mesh->surface->nIDs->datatarray().data();
    if (info::mesh->dimension == 2) {
        std::vector<double> coords; coords.reserve(info::mesh->dimension * size);
        for (size_t n = 0; n < size; ++n) {
            coords.push_back(info::mesh->surface->coordinates->datatarray()[n].x);
            coords.push_back(info::mesh->surface->coordinates->datatarray()[n].y);
        }
        participant.setMeshVertices(info::ecf->coupling.dummy.mesh, coords, precice::span(ids, ids + size));
    } else {
        double *coords = &info::mesh->surface->coordinates->datatarray().data()->x;
        participant.setMeshVertices(info::ecf->coupling.dummy.mesh, precice::span(coords, coords + info::mesh->dimension * size), precice::span(ids, ids + size));
    }

    std::vector<double> readData(size * info::mesh->dimension);
    std::vector<double> writeData(size * info::mesh->dimension);

    for (size_t i = info::mesh->dimension - 1; i < size * info::mesh->dimension; i += info::mesh->dimension) {
        writeData[i] = 1;
    }
    if (participant.requiresInitialData()) {
        participant.writeData(info::ecf->coupling.dummy.mesh, info::ecf->coupling.dummy.data_out, precice::span(ids, ids + size), writeData);
    }
    participant.initialize();

    while (participant.isCouplingOngoing()) {
        if (participant.requiresWritingCheckpoint()) {
            std::cout << "DUMMY: Writing iteration checkpoint\n";
        }

        double dt = participant.getMaxTimeStepSize();
        participant.readData(info::ecf->coupling.dummy.mesh, info::ecf->coupling.dummy.data_in, precice::span(ids, ids + size), dt, readData);

        for (size_t i = info::mesh->dimension - 1; i < size * info::mesh->dimension; i += info::mesh->dimension) {
            writeData[i] = readData[i];
        }

        participant.writeData(info::ecf->coupling.dummy.mesh, info::ecf->coupling.dummy.data_out, precice::span(ids, ids + size), writeData);

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
