
#include "w.precice.h"

#include "analysis/assembler/structuralmechanics.h"
#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/surfacestore.h"

#include <iostream>

#ifdef __clang__
#include <experimental/filesystem>
#else
#include <filesystem>
#endif

#ifdef HAVE_PRECICE
#include <precice/precice.hpp>

namespace espreso {
struct PreciceWrapper {
    PreciceWrapper(const CouplingConfiguration &configuration)
    : configuration(configuration),
      precice(configuration.solver, std::filesystem::path(info::ecf->ecffile).parent_path().append(configuration.configuration).c_str(), info::mpi::rank, info::mpi::size)
    {
        surface = info::mesh->bregion("SURFACE");
    }

    const CouplingConfiguration &configuration;
    precice::Participant precice;

    BoundaryRegionStore *surface;

    struct Mesh {
        std::map<std::string, BoundaryElementData* > data;
        std::vector<precice::VertexID> ids;
        std::vector<double> buffer; // used for exchange other data via PreCICE
    };

    std::map<std::string, Mesh> mesh;
};

static void setMesh2D(PreciceWrapper &w)
{
    Point min(1e20, 1e20, 0), max(-1e20, -1e20, 0);

    if (w.configuration.mesh.size()) {
        size_t size = info::mesh->surface->nIDs->datatarray().size();

        std::vector<double> coords;
        coords.reserve(info::mesh->dimension * size);
        for (size_t n = 0; n < size; ++n) {
            info::mesh->surface->coordinates->datatarray()[n].minmax(min, max);
            coords.push_back(info::mesh->surface->coordinates->datatarray()[n].x);
            coords.push_back(info::mesh->surface->coordinates->datatarray()[n].y);
        }
        w.mesh[w.configuration.mesh].ids.resize(size);
        w.precice.setMeshVertices(w.configuration.mesh, coords, w.mesh[w.configuration.mesh].ids);

        std::vector<precice::VertexID> edges;
        for (auto e = info::mesh->surface->enodes->cbegin(); e != info::mesh->surface->enodes->cend(); ++e) {
            for (auto n = e->begin(); n != e->end(); ++n) {
                edges.push_back(*n);
            }
        }
         w.precice.setMeshEdges(w.configuration.mesh, edges);
    }

     if (w.configuration.centers.size()) {
         size_t size = info::mesh->surface->enodes->structures();

         std::vector<double> coords;
         coords.reserve(info::mesh->dimension * size);
         for (auto e = info::mesh->surface->enodes->cbegin(); e != info::mesh->surface->enodes->cend(); ++e) {
             Point c;
             for (auto n = e->begin(); n != e->end(); ++n) {
                 c += info::mesh->surface->coordinates->datatarray()[*n];
             }
             c /= e->size();
             c.minmax(min, max);
             coords.push_back(c.x);
             coords.push_back(c.y);
         }
         w.mesh[w.configuration.centers].ids.resize(size);
         w.precice.setMeshVertices(w.configuration.centers, coords, w.mesh[w.configuration.centers].ids);
     }


     for (auto data = w.configuration.exchange.cbegin(); data != w.configuration.exchange.cend(); ++data) {
         if (data->second.direct) {
             printf("add %s\n", data->first.c_str());
             w.precice.setMeshAccessRegion(data->first, std::vector<double>{ min.x, max.x, min.y, max.y });
         }
     }
}

static void getMesh2D(PreciceWrapper &w)
{
    PreciceWrapper::Mesh *nodes = nullptr;
    PreciceWrapper::Mesh *centers = nullptr;

    for (auto data = w.configuration.exchange.cbegin(); data != w.configuration.exchange.cend(); ++data) {
        if (data->second.direct) {
            w.mesh[data->first].ids.resize(w.precice.getMeshVertexSize(data->first));
            w.mesh[data->first].buffer.resize(w.mesh[data->first].ids.size() * 2);
            w.precice.getMeshVertexIDsAndCoordinates(data->first, w.mesh[data->first].ids, w.mesh[data->first].buffer);

            if (data->second.centers) {
                if (centers != nullptr) eslog::error("only one center mesh can be set.\n");
                centers = &w.mesh[data->first];
            } else {
                if (nodes != nullptr) eslog::error("only one node mesh can be set.\n");
                nodes = &w.mesh[data->first];
            }
        }
    }


}

static void setMesh3D(PreciceWrapper &w)
{

}

static void getMesh3D(PreciceWrapper &w)
{

}

static void readData(PreciceWrapper &w, double dt)
{
    auto read = [&] (const std::string &mesh, const std::string &data) {
        w.precice.readData(mesh, data, w.mesh[mesh].ids, dt, w.mesh[mesh].data[data]->data);
    };

    auto readDirect = [&] (const std::string &mesh, const std::string &data) {
        w.precice.readData(mesh, data, w.mesh[mesh].ids, dt, w.mesh[mesh].buffer);
        // mortar mapping
    };

    for (auto data = w.configuration.exchange.cbegin(); data != w.configuration.exchange.cend(); ++data) {
        if (0 && data->second.direct) { // TODO: mortar mapping
            if (data->second.read.force)    { readDirect(data->first, "Force"); }
            if (data->second.read.pressure) { readDirect(data->first, "Pressure"); }
            if (data->second.read.stress)   { readDirect(data->first, "Stress"); }
        } else {
            if (data->second.read.force)    { read(data->first, "Force"); }
            if (data->second.read.pressure) { read(data->first, "Pressure"); }
            if (data->second.read.stress)   { read(data->first, "Stress"); }
        }
    }
}

static void writeData(PreciceWrapper &w)
{
    auto write = [&] (const std::string &mesh, const std::string &data) {
        w.precice.writeData(mesh, data, w.mesh[mesh].ids, w.mesh[mesh].data[data]->data);
    };

    for (auto data = w.configuration.exchange.cbegin(); data != w.configuration.exchange.cend(); ++data) {
        if (data->second.write.displacement) { write(data->first, "Displacement");}
        if (data->second.write.velocity)     { write(data->first, "Velocity"); }
    }
}

}

#endif

using namespace espreso;

Precice::Precice(const CouplingConfiguration &configuration)
: wrapper(nullptr)
{
    if (!configuration.isactive()) { return; }

#ifdef HAVE_PRECICE
    wrapper = new PreciceWrapper(configuration);

    switch (info::mesh->dimension) {
    case 2: setMesh2D(*wrapper); break;
    case 3: setMesh3D(*wrapper); break;
    }

    wrapper->precice.initialize();

    switch (info::mesh->dimension) {
    case 2: getMesh2D(*wrapper); break;
    case 3: getMesh3D(*wrapper); break;
    }

    auto resize = [&] (const std::string &mesh, const std::string &data, int dim, bool center) {
        BoundaryElementData::Type type = center ? BoundaryElementData::Type::ELEMENTS : BoundaryElementData::Type::NODES;
        NamedData::DataType datatype = dim == 1 ? NamedData::DataType::SCALAR : NamedData::DataType::VECTOR;
        wrapper->mesh[mesh].data[data] = wrapper->surface->appendData(dim, type, datatype, data);
    };

    for (auto data = configuration.exchange.cbegin(); data != configuration.exchange.cend(); ++data) {
        if (data->second.read.force)         { resize(data->first, "Force"       , info::mesh->dimension, data->second.centers); }
        if (data->second.read.pressure)      { resize(data->first, "Pressure"    ,                     1, data->second.centers); }
        if (data->second.read.stress)        { resize(data->first, "Stress"      , info::mesh->dimension, data->second.centers); }
        if (data->second.write.displacement) { resize(data->first, "Displacement", info::mesh->dimension, data->second.centers); }
        if (data->second.write.velocity)     { resize(data->first, "Velocity"    , info::mesh->dimension, data->second.centers); }
    }
#endif
}

Precice::~Precice()
{
    if (wrapper == nullptr) { return; }
#ifdef HAVE_PRECICE
    delete wrapper;
#endif
}

double Precice::timeStep(double dt)
{
    if (wrapper == nullptr) { return dt; }
#ifdef HAVE_PRECICE
    return std::min(dt, wrapper->precice.getMaxTimeStepSize());
#endif
}

bool Precice::requiresWritingCheckpoint()
{
    if (wrapper == nullptr) { return false; }
#ifdef HAVE_PRECICE
    return wrapper->precice.requiresWritingCheckpoint();
#endif
}

bool Precice::requiresReadingCheckpoint()
{
    if (wrapper == nullptr) { return false; }
#ifdef HAVE_PRECICE
    return wrapper->precice.requiresReadingCheckpoint();
#endif
}

void Precice::read(double dt)
{
    if (wrapper == nullptr) { return; }
#ifdef HAVE_PRECICE
    readData(*wrapper, dt);
#endif
}

void Precice::write()
{
    if (wrapper == nullptr) { return; }
#ifdef HAVE_PRECICE
    writeData(*wrapper);
#endif
}

void Precice::advance(double dt)
{
    if (wrapper == nullptr) { return; }
#ifdef HAVE_PRECICE
    wrapper->precice.advance(dt);
#endif
}

void Precice::dummy()
{
    if (wrapper == nullptr) { return; }
#ifdef HAVE_PRECICE
    auto _read = [&] () {
//        wrapper.precice.writeData(mesh, data, w.mesh[mesh].ids, w.mesh[mesh].data[data]->data);
    };

    while (wrapper->precice.isCouplingOngoing()) {
        if (wrapper->precice.requiresWritingCheckpoint()) {
            // checkpoint
        }
        double dt = wrapper->precice.getMaxTimeStepSize();

//        for (auto data = wrapper->configuration.exchange.cbegin(); data != wrapper->configuration.exchange.cend(); ++data) {
//            if (data->second.read.force)         { resize(data->first, "Force"       , info::mesh->dimension, data->second.centers); }
//            if (data->second.read.pressure)      { resize(data->first, "Pressure"    ,                     1, data->second.centers); }
//            if (data->second.read.stress)        { resize(data->first, "Stress"      , info::mesh->dimension, data->second.centers); }
//            if (data->second.write.displacement) { resize(data->first, "Displacement", info::mesh->dimension, data->second.centers); }
//            if (data->second.write.velocity)     { resize(data->first, "Velocity"    , info::mesh->dimension, data->second.centers); }
//        }
        wrapper->precice.advance(dt);

        if (wrapper->precice.requiresReadingCheckpoint()) {
            // checkpoint
        } else {
            // advance
        }
    }
#endif
}
