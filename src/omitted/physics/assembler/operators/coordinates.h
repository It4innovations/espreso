
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_

#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"

#include "gausspoints.h"
#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"

namespace espreso {

struct CoordinatesToElementNodes: public Operator {
    CoordinatesToElementNodes(serializededata<esint, esint>::const_iterator procNodes, ParameterData &ncoordinates, int interval)
    : Operator(interval, ncoordinates.isconst[interval], ncoordinates.update[interval]),
      procNodes(procNodes),
      ncoordinates(ncoordinates, interval, ndim)
    {

    }

    serializededata<esint, esint>::const_iterator procNodes;
    OutputParameterIterator ncoordinates;

    void operator++()
    {
        ++procNodes;
    }

    CoordinatesToElementNodes& operator+=(const int rhs)
    {
        procNodes += rhs;
        return *this;
    }
};

struct Coordinates2DToElementNodes: CoordinatesToElementNodes {
    using CoordinatesToElementNodes::CoordinatesToElementNodes;

    void operator()()
    {
        for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
            ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
            ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
        }
    }
};

struct Coordinates2DToElementNodesSimd: CoordinatesToElementNodes {

    using CoordinatesToElementNodes::CoordinatesToElementNodes;
    void operator()()
    {
        serializededata<esint, esint>::const_iterator tmpNodeIter(procNodes);
        OutputParameterIterator tmpNcoordinates(ncoordinates);

        for(size_t simdLane = 0; simdLane < SIMD::size && tmpNodeIter->begin() != tmpNodeIter->end(); ++simdLane, ++tmpNodeIter)
        {
            for (auto n = tmpNodeIter->begin(); n != tmpNodeIter->end(); ++n, ++ncoordinates) {
                auto nodeOfElement  = (n - tmpNodeIter->begin());
                tmpNcoordinates[(nodeOfElement * 2 + 0) * SIMD::size + simdLane] = info::mesh->nodes->coordinates->datatarray()[*n].x;
                tmpNcoordinates[(nodeOfElement * 2 + 1) * SIMD::size + simdLane] = info::mesh->nodes->coordinates->datatarray()[*n].y;
            }
        }
    }
};

struct Coordinates3DToElementNodes: CoordinatesToElementNodes {
    using CoordinatesToElementNodes::CoordinatesToElementNodes;

    void operator()()
    {
        for (auto n = procNodes->begin(); n != procNodes->end(); ++n, ++ncoordinates) {
            ncoordinates[0] = info::mesh->nodes->coordinates->datatarray()[*n].x;
            ncoordinates[1] = info::mesh->nodes->coordinates->datatarray()[*n].y;
            ncoordinates[2] = info::mesh->nodes->coordinates->datatarray()[*n].z;
        }
    }
};

struct Coordinates3DToElementNodesSimd: CoordinatesToElementNodes {
    using CoordinatesToElementNodes::CoordinatesToElementNodes;

    void operator()()
    {
        serializededata<esint, esint>::const_iterator tmpNodeIter(procNodes);
        OutputParameterIterator tmpNcoordinates(ncoordinates);

        for(size_t simdLane = 0; simdLane < SIMD::size && tmpNodeIter->begin() != tmpNodeIter->end(); ++simdLane, ++tmpNodeIter)
        {
            for (auto n = tmpNodeIter->begin(); n != tmpNodeIter->end(); ++n, ++ncoordinates) {
                auto nodeOfElement  = (n - tmpNodeIter->begin());
                tmpNcoordinates[(nodeOfElement * 3 + 0) * SIMD::size + simdLane] = info::mesh->nodes->coordinates->datatarray()[*n].x;
                tmpNcoordinates[(nodeOfElement * 3 + 1) * SIMD::size + simdLane] = info::mesh->nodes->coordinates->datatarray()[*n].y;
                tmpNcoordinates[(nodeOfElement * 3 + 2) * SIMD::size + simdLane] = info::mesh->nodes->coordinates->datatarray()[*n].z;
            }
        }
    }
};

struct ElementCoordinates: public ElementOperatorBuilder {
    HeatTransferModuleOpt &kernel;

    ElementCoordinates(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("ELEMENTS COORDINATES"), kernel(kernel)
    {

    }

    bool build(HeatTransferModuleOpt &kernel) override
    {
        kernel.coords.node.addInput(info::mesh->nodes->coordinates);
        kernel.coords.gp.addInput(kernel.coords.node);
        kernel.coords.node.resize();
        kernel.coords.gp.resize();

        kernel.addParameter(kernel.coords.node);
        kernel.addParameter(kernel.coords.gp);

        return true;
    }

    void apply(int interval)
    {
        auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
        if (info::mesh->dimension == 2) {
            iterate_elements(Coordinates2DToElementNodes(procNodes, kernel.coords.node, interval));
            iterate_elements_gps<HeatTransferModuleOpt::NGP>(FromNodesToGaussPoints<2>(kernel.integration.N, kernel.coords.node, kernel.coords.gp, interval));
        }
        if (info::mesh->dimension == 3) {
            iterate_elements(Coordinates3DToElementNodes(procNodes, kernel.coords.node, interval));
            iterate_elements_gps<HeatTransferModuleOpt::NGP>(FromNodesToGaussPoints<3>(kernel.integration.N, kernel.coords.node, kernel.coords.gp, interval));
        }
    }
};

struct ElementCoordinatesSimd: public ElementOperatorBuilder {
    HeatTransferModuleOpt &kernel;

    ElementCoordinatesSimd(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("ELEMENTS COORDINATES"), kernel(kernel)
    {

    }

    bool build(HeatTransferModuleOpt &kernel) override
    {
        kernel.coordsSimd.node.addInput(info::mesh->nodes->coordinates);
        kernel.coordsSimd.gp.addInput(kernel.coordsSimd.node);
        kernel.coordsSimd.node.resizeAligned(SIMD::size * sizeof(double));
        kernel.coordsSimd.gp.resizeAligned(SIMD::size * sizeof(double));

        kernel.addParameter(kernel.coordsSimd.node);
        kernel.addParameter(kernel.coordsSimd.gp);

        return true;
    }

    void apply(int interval)
    {
        auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
        if (info::mesh->dimension == 2) {
            iterate_elements_simd(Coordinates2DToElementNodesSimd(procNodes, kernel.coordsSimd.node, interval));
            iterate_elements_gps_simd<HeatTransferModuleOpt::NGP>(FromNodesToGaussPointsSimd<2>(kernel.integrationSimd.N, kernel.coordsSimd.node, kernel.coordsSimd.gp, interval));
        }
        if (info::mesh->dimension == 3) {
            iterate_elements_simd(Coordinates3DToElementNodesSimd(procNodes, kernel.coordsSimd.node, interval));
            iterate_elements_gps_simd<HeatTransferModuleOpt::NGP>(FromNodesToGaussPointsSimd<3>(kernel.integrationSimd.N, kernel.coordsSimd.node, kernel.coordsSimd.gp, interval));
        }
    }
};

struct BoundaryCoordinates: public BoundaryOperatorBuilder {
    HeatTransferModuleOpt &kernel;

    BoundaryCoordinates(HeatTransferModuleOpt &kernel): BoundaryOperatorBuilder("BOUNDARY COORDINATES"), kernel(kernel)
    {

    }

    bool build(HeatTransferModuleOpt &kernel) override
    {
        for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
            kernel.coords.boundary.node.regions[r].addInput(info::mesh->nodes->coordinates);
            kernel.coords.boundary.gp.regions[r].addInput(kernel.coords.boundary.node.regions[r]);
            kernel.coords.boundary.node.regions[r].isset = true;
            kernel.coords.boundary.gp.regions[r].isset = true;
            kernel.coords.boundary.node.regions[r].resize();
            kernel.coords.boundary.gp.regions[r].resize();

            kernel.addParameter(kernel.coords.boundary.node.regions[r]);
            kernel.addParameter(kernel.coords.boundary.gp.regions[r]);
        }
        return true;
    }

    void apply(int region, int interval)
    {
        auto procNodes = info::mesh->boundaryRegions[region]->nodes->cbegin() + info::mesh->boundaryRegions[region]->eintervals[interval].begin;
        if (info::mesh->dimension == 2) {
            iterate_boundary(Coordinates2DToElementNodes(procNodes, kernel.coords.boundary.node.regions[region], interval), region);
            iterate_boundary_gps<HeatTransferModuleOpt::NGP>(FromNodesToGaussPoints<2>(kernel.integration.boundary.N.regions[region], kernel.coords.boundary.node.regions[region], kernel.coords.boundary.gp.regions[region], interval), region);
        }
        if (info::mesh->dimension == 3) {
            iterate_boundary(Coordinates3DToElementNodes(procNodes, kernel.coords.boundary.node.regions[region], interval), region);
            iterate_boundary_gps<HeatTransferModuleOpt::NGP>(FromNodesToGaussPoints<3>(kernel.integration.boundary.N.regions[region], kernel.coords.boundary.node.regions[region], kernel.coords.boundary.gp.regions[region], interval), region);
        }
    }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_COORDINATES_H_ */
