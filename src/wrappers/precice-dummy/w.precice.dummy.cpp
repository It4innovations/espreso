
#include "w.precice.dummy.h"
#include "basis/containers/serializededata.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/surfacestore.h"

#include <iostream>
#include <vector>

#ifdef HAVE_PRECICE
#include <precice/precice.hpp>
#endif

using namespace espreso;

void PreciceDummy::run()
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

    auto _read  = [&] (const std::string &name, double dt) { participant.readData (info::ecf->coupling.dummy.mesh, name, precice::span(ids, ids + size), dt, readData); };
    auto _write = [&] (const std::string &name)            { participant.writeData(info::ecf->coupling.dummy.mesh, name, precice::span(ids, ids + size), writeData); };

    for (size_t i = info::mesh->dimension - 1; i < size * info::mesh->dimension; i += info::mesh->dimension) {
        writeData[i] = 1;
    }
    if (participant.requiresInitialData()) {
        if (info::ecf->coupling.dummy.data_out.displacement) { _write("Displacement"); }
        if (info::ecf->coupling.dummy.data_out.velocity)     { _write("Velocity"); }
    }
    participant.initialize();

    while (participant.isCouplingOngoing()) {
        if (participant.requiresWritingCheckpoint()) {
            std::cout << "DUMMY: Writing iteration checkpoint\n";
        }

        double dt = participant.getMaxTimeStepSize();
        if (info::ecf->coupling.dummy.data_in.force)    { _read("Force"   , dt); }
        if (info::ecf->coupling.dummy.data_in.pressure) { _read("Pressure", dt); }
        if (info::ecf->coupling.dummy.data_in.stress)   { _read("Stress"  , dt); }

        for (size_t i = info::mesh->dimension - 1; i < size * info::mesh->dimension; i += info::mesh->dimension) {
            writeData[i] = readData[i];
        }

        if (info::ecf->coupling.dummy.data_out.displacement) { _write("Displacement"); }
        if (info::ecf->coupling.dummy.data_out.velocity)     { _write("Velocity"); }

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
