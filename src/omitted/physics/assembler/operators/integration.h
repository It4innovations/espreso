
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_

#include "physics/assembler/operator.h"
#include "physics/assembler/parameter.h"
#include "physics/assembler/math.hpp"
#include "math/simd/simd.h"

namespace espreso {

struct ElementJacobian: public Operator {
    ElementJacobian(
            const ParameterData &coordinates,
            const ParameterData &dN,
            ParameterData &inversion,
            ParameterData &det,
            ParameterData &dND,
            int interval)
    : Operator(interval, false, inversion.update[interval] || det.update[interval] || dND.update[interval]),
      coords(coordinates, interval),
      dN(dN, interval, 0),
      inv(inversion, interval),
      det(det, interval),
      dND(dND, interval)
    {

    }

    InputParameterIterator coords, dN;
    OutputParameterIterator inv, det, dND;

    void operator++()
    {
        ++coords;
        ++inv; ++det; ++dND;
    }

    ElementJacobian& operator+=(const size_t rhs)
    {
        coords +=rhs;
        inv += rhs; det += rhs; dND += rhs;
        return *this;
    }
};

struct ElementJacobian2D: public ElementJacobian {
    using ElementJacobian::ElementJacobian;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double jacobian[4] = { 0, 0, 0, 0 };

        for (int n = 0; n < nodes; ++n) {
            jacobian[0] += dN[2 * gpindex * nodes + n + 0 * nodes] * coords[2 * n + 0];
            jacobian[1] += dN[2 * gpindex * nodes + n + 0 * nodes] * coords[2 * n + 1];
            jacobian[2] += dN[2 * gpindex * nodes + n + 1 * nodes] * coords[2 * n + 0];
            jacobian[3] += dN[2 * gpindex * nodes + n + 1 * nodes] * coords[2 * n + 1];
        }

        det[gpindex] = jacobian[0] * jacobian[3] - jacobian[1] * jacobian[2];
        double detJx = 1 / det[gpindex];
        inv[4 * gpindex + 0] =   detJx * jacobian[3];
        inv[4 * gpindex + 1] = - detJx * jacobian[1];
        inv[4 * gpindex + 2] = - detJx * jacobian[2];
        inv[4 * gpindex + 3] =   detJx * jacobian[0];

        M22M2N<nodes>(inv.data + 4 * gpindex, dN.data + 2 * gpindex * nodes, dND.data + 2 * gpindex * nodes);
    }
};

struct ElementJacobian2DSimd: public ElementJacobian {
    using ElementJacobian::ElementJacobian;

    template<size_t nodes, size_t gps>
    void operator()(size_t gpindex)
    {
        const double * __restrict__ pCoords = coords.data;
        const double * __restrict__ pDN     = dN.data;

        double * __restrict__ pDet = det.data;
        double * __restrict__ pInv = inv.data;
        double * __restrict__ pDND = dND.data;

        SIMD jacobian0, jacobian1, jacobian2, jacobian3;
        SIMD dNX, dNY, coordsX, coordsY;
        
        #pragma unroll(nodes)
        for (size_t n = 0; n < nodes; ++n) {

            coordsX = load(&pCoords[(n * 2 + 0) * SIMD::size]);
            coordsY = load(&pCoords[(n * 2 + 1) * SIMD::size]);
            dNX = load(&pDN[(2 * gpindex * nodes + n + 0 * nodes) * SIMD::size]);
            dNY = load(&pDN[(2 * gpindex * nodes + n + 1 * nodes) * SIMD::size]);

            jacobian0 = jacobian0 + dNX * coordsX;
            jacobian1 = jacobian1 + dNX * coordsY;
            jacobian2 = jacobian2 + dNY * coordsX;
            jacobian3 = jacobian3 + dNY * coordsY;
        }

        SIMD determinant = jacobian0 * jacobian3 - jacobian1 * jacobian2;
        store(&pDet[gpindex * SIMD::size], determinant);

        SIMD detJx = ones() / determinant;
        store(&pInv[(4 * gpindex + 0) * SIMD::size],  detJx * jacobian3);
        store(&pInv[(4 * gpindex + 1) * SIMD::size], -detJx * jacobian1);
        store(&pInv[(4 * gpindex + 2) * SIMD::size], -detJx * jacobian2);
        store(&pInv[(4 * gpindex + 3) * SIMD::size],  detJx * jacobian0);

        M22M2NSimd<nodes>(&pInv[4 * SIMD::size * gpindex], &pDN[2 * SIMD::size * nodes * gpindex], &pDND[2 * SIMD::size * nodes * gpindex]);
    }
};

struct ElementJacobian3D: public ElementJacobian {
    using ElementJacobian::ElementJacobian;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double jacobian[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        for (int n = 0; n < nodes; ++n) {
            jacobian[0] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 0];
            jacobian[1] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 1];
            jacobian[2] += dN[3 * gpindex * nodes + n + 0 * nodes] * coords[3 * n + 2];
            jacobian[3] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 0];
            jacobian[4] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 1];
            jacobian[5] += dN[3 * gpindex * nodes + n + 1 * nodes] * coords[3 * n + 2];
            jacobian[6] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 0];
            jacobian[7] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 1];
            jacobian[8] += dN[3 * gpindex * nodes + n + 2 * nodes] * coords[3 * n + 2];
        }
        det[gpindex] =
                + jacobian[0] * jacobian[4] * jacobian[8]
                + jacobian[1] * jacobian[5] * jacobian[6]
                + jacobian[2] * jacobian[3] * jacobian[7]
                - jacobian[2] * jacobian[4] * jacobian[6]
                - jacobian[1] * jacobian[3] * jacobian[8]
                - jacobian[0] * jacobian[5] * jacobian[7];

        double detJx = 1 / det[gpindex];
        inv[9 * gpindex + 0] = detJx * ( jacobian[8] * jacobian[4] - jacobian[7] * jacobian[5]);
        inv[9 * gpindex + 1] = detJx * (-jacobian[8] * jacobian[1] + jacobian[7] * jacobian[2]);
        inv[9 * gpindex + 2] = detJx * ( jacobian[5] * jacobian[1] - jacobian[4] * jacobian[2]);
        inv[9 * gpindex + 3] = detJx * (-jacobian[8] * jacobian[3] + jacobian[6] * jacobian[5]);
        inv[9 * gpindex + 4] = detJx * ( jacobian[8] * jacobian[0] - jacobian[6] * jacobian[2]);
        inv[9 * gpindex + 5] = detJx * (-jacobian[5] * jacobian[0] + jacobian[3] * jacobian[2]);
        inv[9 * gpindex + 6] = detJx * ( jacobian[7] * jacobian[3] - jacobian[6] * jacobian[4]);
        inv[9 * gpindex + 7] = detJx * (-jacobian[7] * jacobian[0] + jacobian[6] * jacobian[1]);
        inv[9 * gpindex + 8] = detJx * ( jacobian[4] * jacobian[0] - jacobian[3] * jacobian[1]);

        M33M3N<nodes>(inv.data + 9 * gpindex, dN.data + 3 * gpindex * nodes, dND.data + 3 * gpindex * nodes);
    }
};

struct ElementJacobian3DSimd: public ElementJacobian {
    using ElementJacobian::ElementJacobian;

    template<size_t nodes, size_t gps>
    void operator()(size_t gpindex)
    {
        const double * __restrict__ pCoords = coords.data;
        const double * __restrict__ pDN     = dN.data;

        double * __restrict__ pDet = det.data;
        double * __restrict__ pInv = inv.data;
        double * __restrict__ pDND = dND.data;

        SIMD jacobian0, jacobian1, jacobian2, jacobian3, jacobian4, jacobian5, jacobian6, jacobian7, jacobian8;
        SIMD dNX, dNY, dNZ, coordsX, coordsY, coordsZ;

        #pragma unroll(nodes)
        for (size_t n = 0; n < nodes; ++n) {
            coordsX = load(&pCoords[(n * 3 + 0) * SIMD::size]);
            coordsY = load(&pCoords[(n * 3 + 1) * SIMD::size]);
            coordsZ = load(&pCoords[(n * 3 + 2) * SIMD::size]);
            dNX = load(&pDN[(3 * gpindex * nodes + n + 0 * nodes) * SIMD::size]);
            dNY = load(&pDN[(3 * gpindex * nodes + n + 1 * nodes) * SIMD::size]);
            dNZ = load(&pDN[(3 * gpindex * nodes + n + 2 * nodes) * SIMD::size]);

            jacobian0 = jacobian0 + dNX * coordsX;
            jacobian1 = jacobian1 + dNX * coordsY;
            jacobian2 = jacobian2 + dNX * coordsZ;
            jacobian3 = jacobian3 + dNY * coordsX;
            jacobian4 = jacobian4 + dNY * coordsY;
            jacobian5 = jacobian5 + dNY * coordsZ;
            jacobian6 = jacobian6 + dNZ * coordsX;
            jacobian7 = jacobian7 + dNZ * coordsY;
            jacobian8 = jacobian8 + dNZ * coordsZ;
        }

        SIMD determinant =
                + jacobian0 * jacobian4 * jacobian8
                + jacobian1 * jacobian5 * jacobian6
                + jacobian2 * jacobian3 * jacobian7
                - jacobian2 * jacobian4 * jacobian6
                - jacobian1 * jacobian3 * jacobian8
                - jacobian0 * jacobian5 * jacobian7;

        store(&pDet[gpindex * SIMD::size], determinant);

        SIMD detJx = ones() / determinant;

        store(&pInv[(9 * gpindex + 0) * SIMD::size], detJx * ( jacobian8 * jacobian4 - jacobian7 * jacobian5));
        store(&pInv[(9 * gpindex + 1) * SIMD::size], detJx * (-jacobian8 * jacobian1 + jacobian7 * jacobian2));
        store(&pInv[(9 * gpindex + 2) * SIMD::size], detJx * ( jacobian5 * jacobian1 - jacobian4 * jacobian2));
        store(&pInv[(9 * gpindex + 3) * SIMD::size], detJx * (-jacobian8 * jacobian3 + jacobian6 * jacobian5));
        store(&pInv[(9 * gpindex + 4) * SIMD::size], detJx * ( jacobian8 * jacobian0 - jacobian6 * jacobian2));
        store(&pInv[(9 * gpindex + 5) * SIMD::size], detJx * (-jacobian5 * jacobian0 + jacobian3 * jacobian2));
        store(&pInv[(9 * gpindex + 6) * SIMD::size], detJx * ( jacobian7 * jacobian3 - jacobian6 * jacobian4));
        store(&pInv[(9 * gpindex + 7) * SIMD::size], detJx * (-jacobian7 * jacobian0 + jacobian6 * jacobian1));
        store(&pInv[(9 * gpindex + 8) * SIMD::size], detJx * ( jacobian4 * jacobian0 - jacobian3 * jacobian1));

        M33M3NSimd<nodes>(&pInv[9 * gpindex * SIMD::size], &pDN[3 * gpindex * nodes * SIMD::size], &pDND[3 * gpindex * nodes * SIMD::size]);
    }
};

struct ElementIntegration: public ElementOperatorBuilder {
    HeatTransferModuleOpt &kernel;

    ElementIntegration(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("ELEMENTS INTEGRATION"), kernel(kernel)
    {

    }

    bool build(HeatTransferModuleOpt &kernel) override
    {
        kernel.integration.jacobiInversion.addInput(kernel.coords.node);
        kernel.integration.jacobiInversion.addInput(kernel.integration.dN);
        kernel.integration.jacobiInversion.resize();
        kernel.addParameter(kernel.integration.jacobiInversion);

        kernel.integration.jacobiDeterminant.addInput(kernel.coords.node);
        kernel.integration.jacobiDeterminant.addInput(kernel.integration.dN);
        kernel.integration.jacobiDeterminant.resize();
        kernel.addParameter(kernel.integration.jacobiDeterminant);

        kernel.integration.dND.addInput(kernel.coords.node);
        kernel.integration.dND.addInput(kernel.integration.dN);
        kernel.integration.dND.resize();
        kernel.addParameter(kernel.integration.dND);
        
        return true;
    }

    void apply(int interval)
    {
        if (info::mesh->dimension == 2) {
            iterate_elements_gps<HeatTransferModuleOpt::NGP>(ElementJacobian2D(kernel.coords.node, kernel.integration.dN, kernel.integration.jacobiInversion, kernel.integration.jacobiDeterminant, kernel.integration.dND, interval));
        }
        if (info::mesh->dimension == 3) {
            iterate_elements_gps<HeatTransferModuleOpt::NGP>(ElementJacobian3D(kernel.coords.node, kernel.integration.dN, kernel.integration.jacobiInversion, kernel.integration.jacobiDeterminant, kernel.integration.dND, interval));
        }
    }
};

struct ElementIntegrationSimd: public ElementOperatorBuilder {
    HeatTransferModuleOpt &kernel;

    ElementIntegrationSimd(HeatTransferModuleOpt &kernel): ElementOperatorBuilder("ELEMENTS INTEGRATION"), kernel(kernel)
    {

    }

    bool build(HeatTransferModuleOpt &kernel) override
    {
        kernel.integrationSimd.jacobiInversion.addInput(kernel.coordsSimd.node);
        kernel.integrationSimd.jacobiInversion.addInput(kernel.integration.dN);
        kernel.integrationSimd.jacobiInversion.resizeAligned(SIMD::size * sizeof(double));
        kernel.addParameter(kernel.integrationSimd.jacobiInversion);

        kernel.integrationSimd.jacobiDeterminant.addInput(kernel.coordsSimd.node);
        kernel.integrationSimd.jacobiDeterminant.addInput(kernel.integration.dN);
        kernel.integrationSimd.jacobiDeterminant.resizeAligned(SIMD::size * sizeof(double));
        kernel.addParameter(kernel.integrationSimd.jacobiDeterminant);

        kernel.integrationSimd.dND.addInput(kernel.coordsSimd.node);
        kernel.integrationSimd.dND.addInput(kernel.integration.dN);
        kernel.integrationSimd.dND.resizeAligned(SIMD::size * sizeof(double));
        kernel.addParameter(kernel.integrationSimd.dND);

        return true;
    }

    void apply(int interval)
    {
        if (info::mesh->dimension == 2) {
            iterate_elements_gps_simd<HeatTransferModuleOpt::NGP>(ElementJacobian2DSimd(kernel.coordsSimd.node, kernel.integrationSimd.dN, kernel.integrationSimd.jacobiInversion, kernel.integrationSimd.jacobiDeterminant, kernel.integrationSimd.dND, interval));
        }
        if (info::mesh->dimension == 3) {
            iterate_elements_gps_simd<HeatTransferModuleOpt::NGP>(ElementJacobian3DSimd(kernel.coordsSimd.node, kernel.integrationSimd.dN, kernel.integrationSimd.jacobiInversion, kernel.integrationSimd.jacobiDeterminant, kernel.integrationSimd.dND, interval));
        }
    }
};

struct BoundaryJacobian: public Operator {
    BoundaryJacobian(
            const ParameterData &coordinates,
            const ParameterData &dN,
            ParameterData &jacobian,
            int interval)
    : Operator(interval, jacobian.isconst[interval], jacobian.update[interval]),
      coords(coordinates, interval),
      dN(dN, interval, 0),
      jacobian(jacobian, interval)
    { }

    InputParameterIterator coords, dN;
    OutputParameterIterator jacobian;

    void operator++()
    {
        ++coords;
        ++jacobian;
    }
};

struct BoundaryFaceJacobian: public BoundaryJacobian {
    using BoundaryJacobian::BoundaryJacobian;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double dND[6] = { 0, 0, 0, 0, 0, 0 };
        M2NMN3<nodes>(1, dN.data + 2 * gpindex * nodes, coords.data, dND);
        double x = dND[1] * dND[5] - dND[2] * dND[4];
        double y = dND[2] * dND[3] - dND[0] * dND[5];
        double z = dND[0] * dND[4] - dND[1] * dND[3];
        jacobian.data[gpindex] = std::sqrt(x * x + y * y + z * z);
    }
};

struct BoundaryEdge3DJacobian: public BoundaryJacobian {
    using BoundaryJacobian::BoundaryJacobian;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double dND[3] = { 0, 0, 0 };
        M1NMN3<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
        jacobian.data[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1] + dND[2] * dND[2]);
    }
};

struct BoundaryEdge2DJacobian: public BoundaryJacobian {
    using BoundaryJacobian::BoundaryJacobian;

    template<int nodes, int gps>
    void operator()(int gpindex)
    {
        double dND[2] = { 0, 0 };
        M1NMN2<nodes>(1, dN.data + 1 * gpindex * nodes, coords.data, dND);
        jacobian.data[gpindex] = std::sqrt(dND[0] * dND[0] + dND[1] * dND[1]);
    }
};

struct BoundaryIntegration: public BoundaryOperatorBuilder {
    HeatTransferModuleOpt &kernel;

    BoundaryIntegration(HeatTransferModuleOpt &kernel): BoundaryOperatorBuilder("BOUNDARY INTERGRATION"), kernel(kernel)
    {

    }

    bool build(HeatTransferModuleOpt &kernel) override
    {
        for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
            kernel.integration.boundary.jacobian.regions[r].addInput(kernel.coords.boundary.node.regions[r]);
            kernel.integration.boundary.jacobian.regions[r].addInput(kernel.integration.boundary.dN.regions[r]);
            kernel.integration.boundary.jacobian.regions[r].isset = true;
            kernel.integration.boundary.jacobian.regions[r].resize();
            kernel.addParameter(kernel.integration.boundary.jacobian.regions[r]);
        }
        return true;
    }

    void apply(int region, int interval)
    {
        if (info::mesh->boundaryRegions[region]->dimension == 2) {
            iterate_boundary_gps<HeatTransferModuleOpt::NGP>(BoundaryFaceJacobian(kernel.coords.boundary.node.regions[region], kernel.integration.boundary.dN.regions[region], kernel.integration.boundary.jacobian.regions[region], interval), region);
        }
        if (info::mesh->boundaryRegions[region]->dimension == 1) {
            if (info::mesh->dimension == 3) {
                iterate_boundary_gps<HeatTransferModuleOpt::NGP>(BoundaryEdge3DJacobian(kernel.coords.boundary.node.regions[region], kernel.integration.boundary.dN.regions[region], kernel.integration.boundary.jacobian.regions[region], interval), region);
            }
            if (info::mesh->dimension == 2) {
                iterate_boundary_gps<HeatTransferModuleOpt::NGP>(BoundaryEdge2DJacobian(kernel.coords.boundary.node.regions[region], kernel.integration.boundary.dN.regions[region], kernel.integration.boundary.jacobian.regions[region], interval), region);
            }
        }
    }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_ */
