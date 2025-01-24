
#include "w.precice.h"

#include "analysis/assembler/structuralmechanics.h"
#include "basis/containers/serializededata.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "mesh/store/surfacestore.h"

#include <iostream>
#include <filesystem>

#ifdef HAVE_PRECICE
#include <precice/precice.hpp>

namespace espreso {
struct PreciceData {
    PreciceData()
    : precice(info::ecf->coupling.solver, std::filesystem::path(info::ecf->ecffile).parent_path().append(info::ecf->coupling.configuration).c_str(), info::mpi::rank, info::mpi::size),
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
    if (info::ecf->coupling.isactive()) {
        _data = new PreciceData();

        _data->size = info::mesh->surface->nIDs->datatarray().size();
        _data->ids.resize(_data->size);

        if (info::mesh->dimension == 2) {
            std::vector<double> coords; coords.reserve(info::mesh->dimension * _data->size);
            for (size_t n = 0; n < _data->size; ++n) {
                coords.push_back(info::mesh->surface->coordinates->datatarray()[n].x);
                coords.push_back(info::mesh->surface->coordinates->datatarray()[n].y);
            }
            _data->precice.setMeshVertices(info::ecf->coupling.mesh, coords, _data->ids);
            std::vector<esint> edges;
            for (auto e = info::mesh->surface->enodes->cbegin(); e != info::mesh->surface->enodes->cend(); ++e) {
                for (auto n = e->begin(); n != e->end(); ++n) {
                    edges.push_back(*n);
                }
            }
            _data->precice.setMeshEdges(info::ecf->coupling.mesh, edges);
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

void Precice::_read(double *data, const std::string &name, double dt)
{
#ifdef HAVE_PRECICE
    _data->precice.readData(info::ecf->coupling.mesh, name, _data->ids, dt, _data->data);
    for (size_t n = 0; n < _data->size; ++n) {
        for (int d = 0; d < info::mesh->dimension; ++d) {
            data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d] = _data->data[n * info::mesh->dimension + d];
        }
    }
#endif
}

void Precice::read(double dt)
{
#ifdef HAVE_PRECICE
    if (_data) {
        if (info::ecf->coupling.data_in.force)    { _read(StructuralMechanics::Results::fluidForce->data.data()   , "Force"   , dt); }
        if (info::ecf->coupling.data_in.pressure) { _read(StructuralMechanics::Results::fluidPressure->data.data(), "Pressure", dt); }
        if (info::ecf->coupling.data_in.stress)   { _read(StructuralMechanics::Results::fluidStress->data.data()  , "Stress"  , dt); }
    }
#endif
}

void Precice::_write(double *data, const std::string &name)
{
#ifdef HAVE_PRECICE
    for (size_t n = 0; n < _data->size; ++n) {
        for (int d = 0; d < info::mesh->dimension; ++d) {
            _data->data[n * info::mesh->dimension + d] = data[info::mesh->surface->nIDs->datatarray()[n] * info::mesh->dimension + d];
        }
    }
    _data->precice.writeData(info::ecf->coupling.mesh, name, _data->ids, _data->data);
#endif
}

void Precice::write()
{
#ifdef HAVE_PRECICE
    if (_data) {
        if (info::ecf->coupling.data_out.displacement) { _write(StructuralMechanics::Results::displacement->data.data(), "Displacement"); }
        if (info::ecf->coupling.data_out.velocity)     { _write(StructuralMechanics::Results::velocity->data.data()    , "Velocity"); }
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
