
#include "neutralmesh.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/utils.h"
#include "esinfo/eslog.h"
#include "mesh/element.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"

#include <numeric>

using namespace espreso;

NetgenNeutralMesh::NetgenNeutralMesh(InputFilePack &meshfile)
: _meshfile(meshfile)
{

}

void NetgenNeutralMesh::parse(MeshBuilder &mesh)
{
    profiler::syncstart("netgen_parse");
    profiler::syncparam("size", _meshfile.end - _meshfile.begin);
    int threads = info::env::OMP_NUM_THREADS;
    std::vector<size_t> tdistribution = tarray<size_t>::distribute(threads, _meshfile.end - _meshfile.begin);
    for (int t = 1; t < threads; t++) {
        const char* c = _meshfile.begin + tdistribution[t];
        if (*(c - 1) != '\n') {
            while (c < _meshfile.end && *c++ != '\n');// start at new line
            tdistribution[t] += c - (_meshfile.begin + tdistribution[t]);
        }
    }

    esint nlines = 0, offset;
    std::vector<esint> tlines(threads);
    std::vector<std::vector<esint> > sizes(threads);

    #pragma omp parallel for
    for (int t = 0; t < threads; t++) {
        esint lines = 0, spaces = 0;
        std::vector<esint> tsizes;

        for (auto c = _meshfile.begin + tdistribution[t], cc = c; c < _meshfile.begin + tdistribution[t + 1]; ++c) {
            if (*c == ' ') { ++spaces; }
            if (*c == '\n') {
                if (spaces == 0) {
                    tsizes.push_back(atol(cc));
                }
                ++lines;
                cc = c + 1;
                spaces = 0;
            }
        }

        tlines[t] = lines;
        sizes[t].swap(tsizes);
    }
    profiler::synccheckpoint("scan");

    nlines += tlines[0];
    for (int t = 1; t < threads; t++) {
        nlines += tlines[t];
        sizes[0].insert(sizes[0].end(), sizes[t].begin(), sizes[t].end());
    }
    utils::sizesToOffsets(tlines);

    offset = nlines;
    Communication::exscan(offset);
    Communication::allGatherUnknownSize(sizes[0]);

    profiler::synccheckpoint("synchronize");

    if (sizes[0].size() != 3) {
        eslog::error("Netgen neutral parser: unknown format of data.\n");
    }

    esint pbegin = 1, pend = pbegin + sizes[0][0];
    esint ebegin = pend + 1, eend = ebegin + sizes[0][1];
    esint bbegin = eend + 1, bend = bbegin + sizes[0][2];

    std::vector<std::vector<Point> > tpoints(threads);
    std::vector<std::vector<esint> > telements(threads), tboundary(threads), tid(threads);
    std::vector<std::vector<std::vector<esint> > > tereg(threads), tbreg(threads);

    #pragma omp parallel for
    for (int t = 0; t < threads; t++) {
        std::vector<Point> points;
        std::vector<esint> elements, boundary;
        std::vector<std::vector<esint> > ereg, breg;
        esint line = offset + tlines[t];
        const char *c = _meshfile.begin + tdistribution[t];
        char *cc;

        while (line < pbegin && c < _meshfile.begin + tdistribution[t + 1]) { if (*c++ == '\n') { ++line; } }
        for (; line < pend && c < _meshfile.begin + tdistribution[t + 1]; ++line, c = ++cc) {
            points.push_back(Point{});
            points.back().x = strtod(c, &cc); points.back().y = strtod(cc, &cc); points.back().z = strtod(cc, &cc);
        }
        while (line < ebegin && c < _meshfile.begin + tdistribution[t + 1]) { if (*c++ == '\n') { ++line; } }
        for (; line < eend && c < _meshfile.begin + tdistribution[t + 1]; ++line, c = ++cc) {
            size_t region = strtoul(c, &cc, 10);
            if (ereg.size() <= region) {
                ereg.resize(region + 1);
            }
            ereg[region].push_back(line - ebegin);
            elements.push_back(strtol(cc, &cc, 10) - 1);
            elements.push_back(strtol(cc, &cc, 10) - 1);
            elements.push_back(strtol(cc, &cc, 10) - 1);
            elements.push_back(strtol(cc, &cc, 10) - 1);
            c = cc;
        }
        while (line < bbegin && c < _meshfile.begin + tdistribution[t + 1]) { if (*c++ == '\n') { ++line; } }
        for (; line < bend && c < _meshfile.begin + tdistribution[t + 1]; ++line, c = ++cc) {
            size_t region = strtoul(c, &cc, 10);
            if (breg.size() <= region) {
                breg.resize(region + 1);
            }
            breg[region].push_back(line - ebegin - 1);
            boundary.push_back(strtol(cc, &cc, 10) - 1);
            boundary.push_back(strtol(cc, &cc, 10) - 1);
            boundary.push_back(strtol(cc, &cc, 10) - 1);
            c = cc;
        }

        tpoints[t].swap(points);
        telements[t].swap(elements);
        tboundary[t].swap(boundary);
        tereg[t].swap(ereg);
        tbreg[t].swap(breg);
    }

    esint ntetra = 0, ntria = 0;
    for (int t = 0; t < threads; t++) {
        mesh.coordinates.insert(mesh.coordinates.end(), tpoints[t].begin(), tpoints[t].end());
        mesh.enodes.insert(mesh.enodes.end(), telements[t].begin(), telements[t].end());
        ntetra += telements[t].size() / 4;
    }
    mesh.nIDs.resize(mesh.coordinates.size());
    esint prevpoints = 0;
    if (pbegin < offset) {
        prevpoints = offset - pbegin;
    }
    std::iota(mesh.nIDs.begin(), mesh.nIDs.end(), prevpoints);
    for (int t = 0; t < threads; t++) {
        mesh.enodes.insert(mesh.enodes.end(), tboundary[t].begin(), tboundary[t].end());
        ntria += tboundary[t].size() / 3;
    }
    mesh.etype.resize(ntetra, (int)Element::CODE::TETRA4);
    mesh.etype.resize(ntetra + ntria, (int)Element::CODE::TRIANGLE3);
    mesh.esize.resize(ntetra, 4);
    mesh.esize.resize(ntetra + ntria, 3);
    mesh.eIDs.resize(mesh.esize.size());
    esint prevelements = 0;
    if (ebegin < offset) {
        prevelements = offset - ebegin;
    }
    if (bbegin < offset) {
        --prevelements;
    }
    std::iota(mesh.eIDs.begin(), mesh.eIDs.end(), prevelements);

    std::vector<size_t> nregs= { 0, 0 };
    for (int t = 0; t < threads; t++) {
        nregs[0] = std::max(nregs[0], tereg[t].size());
        nregs[1] = std::max(nregs[1], tbreg[t].size());
    }

    Communication::allReduce(nregs, Communication::OP::MAX);

    for (size_t r = 1; r < nregs[0]; ++r) {
        auto &ids = mesh.eregions["ELEMENTS_" + std::to_string(r)];
        for (int t = 0; t < threads; t++) {
            if (r < tereg[t].size()) {
                ids.insert(ids.end(), tereg[t][r].begin(), tereg[t][r].end());
            }
        }
    }
    for (size_t r = 1; r < nregs[1]; ++r) {
        auto &ids = mesh.eregions["BOUNDARY_" + std::to_string(r)];
        for (int t = 0; t < threads; t++) {
            if (r < tbreg[t].size()) {
                ids.insert(ids.end(), tbreg[t][r].begin(), tbreg[t][r].end());
            }
        }
    }

    profiler::syncend("netgen_parse");
}
