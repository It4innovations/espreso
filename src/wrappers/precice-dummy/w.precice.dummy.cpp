
#include "w.precice.dummy.h"
#include "basis/containers/serializededata.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/mpiinfo.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "output/output.h"

#include <iostream>
#include <vector>

#ifdef __clang__
#include <experimental/filesystem>
#else
#include <filesystem>
#endif

#ifdef HAVE_PRECICE
#include <precice/precice.hpp>
#endif

using namespace espreso;

void PreciceDummy::run()
{
#ifdef HAVE_PRECICE
//    using namespace precice;
//
//    std::cout << "DUMMY: Running solver dummy with preCICE config file \"" << std::filesystem::path(info::ecf->ecffile).parent_path().append(info::ecf->coupling.configuration).c_str() << "\" and participant name \"" << info::ecf->coupling.dummy.solver << "\".\n";
//
//    Participant participant(info::ecf->coupling.dummy.solver, std::filesystem::path(info::ecf->ecffile).parent_path().append(info::ecf->coupling.configuration).c_str(), info::mpi::rank, info::mpi::size);
//
//
//
//    size_t size = info::mesh->surface->nIDs->datatarray().size();
//    std::vector<esint> ids; // we must copy ids since preCICE rewrite them
//    ids.insert(ids.end(), info::mesh->surface->nIDs->datatarray().begin(), info::mesh->surface->nIDs->datatarray().end());
//    std::vector<double> coords; coords.reserve(info::mesh->dimension * size);
//    for (size_t n = 0; n < size; ++n) {
//        for (int d = 0; d < info::mesh->dimension; ++d) {
//            coords.push_back(info::mesh->surface->coordinates->datatarray()[n][d]);
//        }
//    }
//    participant.setMeshVertices(info::ecf->coupling.dummy.mesh, coords, ids);
//
//    std::vector<esint> cids;
//    size_t csize = info::mesh->surface->enodes->structures();
//    if (info::ecf->coupling.dummy.centers.size()) {
//        size_t csize = info::mesh->surface->enodes->structures();
//        cids.resize(csize);
//        std::vector<double> ccoords; ccoords.reserve(info::mesh->dimension * csize);
//        for (auto e = info::mesh->surface->enodes->cbegin(); e != info::mesh->surface->enodes->cend(); ++e) {
//            Point p;
//            for (auto n = e->begin(); n != e->end(); ++n) {
//                p += info::mesh->surface->coordinates->datatarray()[*n];
//            }
//            p /= e->size();
//
//            for (int d = 0; d < info::mesh->dimension; ++d) {
//                ccoords.push_back(p[d]);
//            }
//        }
//        participant.setMeshVertices(info::ecf->coupling.dummy.centers, ccoords, cids);
//    }
//
//    NodeData *force        = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "Force");
//    NodeData *pressure     = info::mesh->nodes->appendData(1                    , NamedData::DataType::SCALAR, "Pressure");
//    NodeData *stress       = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "Stress");
//    NodeData *displacement = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "Displacement");
//    NodeData *velocity     = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "Velocity");
//    info::mesh->output->updateMesh();
//
//    Point pforce   (0, 1e-2, 0);
//    Point ppressure(0, 1e+0, 0);
//    Point pstress  (0, 1e+1, 0);
//    for (size_t n = 0; n < size; ++n) {
//        double scale = info::mesh->surface->coordinates->datatarray()[n].y;
//        pressure->data[info::mesh->surface->nIDs->datatarray()[n]] = scale * ppressure.y;
//        for (int d = 0; d < info::mesh->dimension; ++d) {
//            force->data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d] = scale * pforce[d];
//            stress->data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d] = scale * pstress[d];
//        }
//    }
//
//    std::vector<double> preciceData(size * info::mesh->dimension);
//    std::vector<double> cpreciceData(csize * info::mesh->dimension);
//
//    auto _read  = [&] (NodeData *data, double dt) {
//        printf("READ: %s\n", data->name.c_str());
//        participant.readData(info::ecf->coupling.dummy.mesh, data->name, ids, dt, preciceData);
//        Point sum;
//        for (size_t n = 0; n < size; ++n) {
//            for (int d = 0; d < info::mesh->dimension; ++d) {
//                data->data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d] = preciceData[n * info::mesh->dimension + d];
//                sum[d] += preciceData[n * info::mesh->dimension + d];
//            }
//        }
////        scale = sum.length();
//    };
//    auto _write = [&] (NodeData *data) {
//        for (size_t n = 0; n < size; ++n) {
//            for (int d = 0; d < info::mesh->dimension; ++d) {
//                preciceData[n * info::mesh->dimension + d] = data->data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d];
//            }
//        }
//        printf("WRITE: %s\n", data->name.c_str());
//        participant.writeData(info::ecf->coupling.dummy.mesh, data->name, ids, preciceData);
//    };
//
//    if (participant.requiresInitialData()) {
//        if (info::ecf->coupling.dummy.data_in.force)    { _write(force); }
//        if (info::ecf->coupling.dummy.data_in.pressure) { _write(pressure); }
//        if (info::ecf->coupling.dummy.data_in.stress)   { _write(stress); }
//    }
//    participant.initialize();
//
//    step::Step step;
//    info::mesh->output->updateMonitors(step);
//    step::Time time;
//    time.current = 0;
//    while (participant.isCouplingOngoing()) {
//        if (participant.requiresWritingCheckpoint()) {
//            std::cout << "DUMMY: Writing iteration checkpoint\n";
//        }
//
//        double dt = participant.getMaxTimeStepSize();
//        if (info::ecf->coupling.dummy.data_out.displacement) { _read(displacement, dt); }
//        if (info::ecf->coupling.dummy.data_out.velocity)     { _read(velocity, dt); }
//
//        if (info::ecf->coupling.dummy.data_in.force)    { _write(force); }
//        if (info::ecf->coupling.dummy.data_in.pressure) { _write(pressure); }
//        if (info::ecf->coupling.dummy.data_in.stress)   { _write(stress); }
//
//        participant.advance(dt);
//
//        if (participant.requiresReadingCheckpoint()) {
//            std::cout << "DUMMY: Reading iteration checkpoint\n";
//        } else {
//            std::cout << "DUMMY: Advancing in time\n";
//        }
//        step.substep += 1;
//        time.current += dt;
//        info::mesh->output->updateSolution(step, time);
//    }
//
//    participant.finalize();
//    std::cout << "DUMMY: Closing C++ solver dummy...\n";
#endif
}
