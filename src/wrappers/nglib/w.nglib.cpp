
#ifdef HAVE_NGLIB
#include <cstddef>
namespace nglib{
#include "nglib.h"
}
using namespace nglib;
#endif

#include "w.nglib.h"
#include "basis/utilities/parser.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/mpiinfo.h"
#include "input/meshbuilder.h"
#include "mesh/element.h"

using namespace espreso;

bool NGLib::islinked()
{
#ifdef HAVE_NGLIB
    return true;
#else
    return false;
#endif
}

void NGLib::generate(MeshData &mesh)
{
#ifndef HAVE_NGLIB
    eslog::globalerror("ESPRESO run-time error: cannot call GMSH library (the library is not linked).\n");
#else
    if (info::mpi::rank) {
        // NGLIB is sequential only
        return;
    }

    auto check = [] (Ng_Result result)
    {
        switch (result) {
        case NG_ERROR: eslog::warning("NG ERROR\n"); break;
        case NG_OK: break;
        case NG_SURFACE_INPUT_ERROR: eslog::warning("NG SURFACE INPUT ERROR\n"); break;
        case NG_VOLUME_FAILURE: eslog::warning("NG VOLUME FAILURE\n"); break;
        case NG_STL_INPUT_ERROR: eslog::warning("NG STL INPUT ERROR\n"); break;
        case NG_SURFACE_FAILURE: eslog::warning("NG SURFACE FAILURE\n"); break;
        case NG_FILE_NOT_FOUND: eslog::warning("NG FILE NOT FOUND\n"); break;
        }
        return result != NG_OK;
    };

    const NGLibConfiguration &configuration = info::ecf->input.generation.nglib_options;
    bool isvalid = false;

    Ng_Meshing_Parameters mp;
    mp.uselocalh = configuration.uselocalh;
    mp.maxh = configuration.maxh;
    mp.minh = configuration.minh;
    mp.fineness = configuration.fineness;
    mp.grading = configuration.grading;
    mp.elementsperedge = configuration.elementsperedge;
    mp.elementspercurve = configuration.elementspercurve;
    mp.closeedgeenable = configuration.closeedgeenable;
    mp.closeedgefact = configuration.closeedgefact;
    mp.minedgelenenable = configuration.minedgelenenable;
    mp.minedgelen = configuration.minedgelen;
    mp.second_order = configuration.second_order;
    mp.quad_dominated = configuration.quad_dominated;
    mp.optsurfmeshenable = configuration.optsurfmeshenable;
    mp.optvolmeshenable = configuration.optvolmeshenable;
    mp.optsteps_3d = configuration.optsteps_3d;
    mp.optsteps_2d = configuration.optsteps_2d;
    mp.invert_tets = configuration.invert_tets;
    mp.invert_trigs = configuration.invert_trigs;
    mp.check_overlap = configuration.check_overlap;
    mp.check_overlapping_boundary = configuration.check_overlapping_boundary;

    Ng_Init();
    Ng_Mesh *ng = Ng_NewMesh();
    mp.Transfer_Parameters();

    if (StringCompare::caseSensitiveSuffix(info::ecf->input.path, "stl")) {
        Ng_STL_Geometry *geom = Ng_STL_LoadGeometry(info::ecf->input.path.c_str());
        check(Ng_STL_InitSTLGeometry(geom));
        check(Ng_STL_MakeEdges(geom, ng, &mp));
        check(Ng_STL_GenerateSurfaceMesh(geom, ng, &mp));
        check(Ng_GenerateVolumeMesh(ng, &mp));
        isvalid = true;
    }

#ifdef OCCGEOMETRY
    if (
            StringCompare::caseSensitiveSuffix(info::ecf->input.path, "stp") ||
            StringCompare::caseSensitiveSuffix(info::ecf->input.path, "step")) {

        Ng_OCC_Geometry *geom = Ng_OCC_Load_STEP(info::ecf->input.path.c_str());
        check(Ng_OCC_SetLocalMeshSize(geom, ng, &mp));
        check(Ng_OCC_GenerateEdgeMesh(geom, ng, &mp));
        check(Ng_OCC_GenerateSurfaceMesh(geom, ng, &mp));
        check(Ng_GenerateVolumeMesh(ng, &mp));
        isvalid = true;
    }

    if (StringCompare::caseSensitiveSuffix(info::ecf->input.path, "igs")) {
        Ng_OCC_Geometry *geom = Ng_OCC_Load_IGES(info::ecf->input.path.c_str());
        check(Ng_OCC_SetLocalMeshSize(geom, ng, &mp));
        check(Ng_OCC_GenerateEdgeMesh(geom, ng, &mp));
        check(Ng_OCC_GenerateSurfaceMesh(geom, ng, &mp));
        check(Ng_GenerateVolumeMesh(ng, &mp));
        isvalid = true;
    }
#endif

    if (!isvalid) {
        eslog::error("MESIO run-time error: NGLib does not support provided format.\n");
    }

    { // insert nodes
        int points = Ng_GetNP(ng);
        double point[3];
        mesh.nIDs.reserve(points);
        mesh.coordinates.reserve(points);
        for (int p = 0; p < points; ++p) {
            Ng_GetPoint(ng, p + 1, point);
            mesh.nIDs.push_back(p);
            mesh.coordinates.push_back(Point(point[0], point[1], point[2]));
        }
    }
    mesh.eIDs.reserve(Ng_GetNE(ng) + Ng_GetNSE(ng));
    mesh.esize.reserve(Ng_GetNE(ng) + Ng_GetNSE(ng));
    mesh.etype.reserve(Ng_GetNE(ng) + Ng_GetNSE(ng));
    mesh.enodes.reserve(Ng_GetNE(ng) * NG_VOLUME_ELEMENT_MAXPOINTS + Ng_GetNSE(ng) * NG_SURFACE_ELEMENT_MAXPOINTS);
    { // insert elements
        int elements = Ng_GetNE(ng);
        int element[NG_VOLUME_ELEMENT_MAXPOINTS];
        for (int e = 0; e < elements; e++) {
            mesh.eIDs.push_back(e);
            switch (Ng_GetVolumeElement(ng, e + 1, element)) {
            case NG_TET:
                mesh.esize.push_back(4);
                mesh.etype.push_back((int)Element::CODE::TETRA4);
                break;
            case NG_PYRAMID:
                mesh.esize.push_back(5);
                mesh.etype.push_back((int)Element::CODE::PYRAMID5);
                break;
            case NG_PRISM:
                mesh.esize.push_back(6);
                mesh.etype.push_back((int)Element::CODE::PRISMA6);
                break;
            case NG_TET10:
                mesh.esize.push_back(10);
                mesh.etype.push_back((int)Element::CODE::TETRA10);
                break;
            }
            for (int i = 0; i < mesh.esize.back(); ++i) {
                mesh.enodes.push_back(element[i] - 1);
            }
        }
        mesh.eregions["VOLUME"] = mesh.eIDs;
    }

    { // insert surface
        int offset = mesh.eIDs.size();
        int elements = Ng_GetNSE(ng);
        int element[NG_SURFACE_ELEMENT_MAXPOINTS];
        for (int e = 0; e < elements; e++) {
            mesh.eIDs.push_back(e + offset);
            switch (Ng_GetSurfaceElement(ng, e + 1, element)) {
            case NG_TRIG:
                mesh.esize.push_back(3);
                mesh.etype.push_back((int)Element::CODE::TRIANGLE3);
                break;
            case NG_QUAD:
                mesh.esize.push_back(4);
                mesh.etype.push_back((int)Element::CODE::SQUARE4);
                break;
            case NG_TRIG6:
                mesh.esize.push_back(6);
                mesh.etype.push_back((int)Element::CODE::TRIANGLE6);
                break;
            case NG_QUAD6:
                eslog::error("MESIO run-time error: NGlib generated unsupported boundary element.\n");
                break;
            case NG_QUAD8:
                mesh.esize.push_back(8);
                mesh.etype.push_back((int)Element::CODE::SQUARE8);
                break;
            }
            for (int i = 0; i < mesh.esize.back(); ++i) {
                mesh.enodes.push_back(element[i] - 1);
            }
        }
        mesh.eregions["SURFACE"].assign(mesh.eIDs.begin() + offset, mesh.eIDs.end());
    }
    Ng_Exit ();
#endif
}
