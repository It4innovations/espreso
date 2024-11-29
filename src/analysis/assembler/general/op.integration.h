
#ifndef SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_
#define SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_

#include "subkernel.h"
#include "math.h"

namespace espreso {

struct Integration: SubKernel {
    Integration()
    {
        isconst = false;
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
    }
};


template <size_t nodes, size_t ndim, size_t edim> struct IntegrationKernel;

template <size_t nodes>
struct IntegrationKernel<nodes, 2, 2>: Integration {
    IntegrationKernel(const Integration &base): Integration(base) {}

    template <typename Element>
    void coords(Element &element, size_t gp)
    {
        SIMD J[4];

        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0];
            SIMD coordsY = element.coords.node[n][1];
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);

            J[0] = J[0] + dNX * coordsX;
            J[2] = J[2] + dNX * coordsY;
            J[1] = J[1] + dNY * coordsX;
            J[3] = J[3] + dNY * coordsY;
        }

        inv22(J, element.det, element.invJ);

        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);
            element.dND[n * 2 + 0] = element.invJ[0] * dNX + element.invJ[2] * dNY;
            element.dND[n * 2 + 1] = element.invJ[1] * dNX + element.invJ[3] * dNY;
        }
    }

    template <typename Element>
    void displacement(Element &element, size_t gp)
    {
        SIMD J[4];

        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);

            J[0] = J[0] + dNX * coordsX;
            J[2] = J[2] + dNX * coordsY;
            J[1] = J[1] + dNY * coordsX;
            J[3] = J[3] + dNY * coordsY;
        }

        inv22(J, element.det, element.invJ);

        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);
            element.dND[n * 2 + 0] = element.invJ[0] * dNX + element.invJ[2] * dNY;
            element.dND[n * 2 + 1] = element.invJ[1] * dNX + element.invJ[3] * dNY;
        }
    }

    template <typename Element>
    void displacementInNodes(Element &element, size_t node)
    {
        SIMD J[4];

        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD dNX = load1(element.dNN[node][n][0]);
            SIMD dNY = load1(element.dNN[node][n][1]);

            J[0] = J[0] + dNX * coordsX;
            J[2] = J[2] + dNX * coordsY;
            J[1] = J[1] + dNY * coordsX;
            J[3] = J[3] + dNY * coordsY;
        }

        inv22(J, element.det, element.invJ);

        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(element.dNN[node][n][0]);
            SIMD dNY = load1(element.dNN[node][n][1]);
            element.dND[n * 2 + 0] = element.invJ[0] * dNX + element.invJ[2] * dNY;
            element.dND[n * 2 + 1] = element.invJ[1] * dNX + element.invJ[3] * dNY;
        }
    }
};

template <size_t nodes>
struct IntegrationKernel<nodes, 3, 3>: Integration {
    IntegrationKernel(const Integration &base): Integration(base) {}

    template <typename Element>
    void coords(Element &element, size_t gp)
    {
        SIMD J[9];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0];
            SIMD coordsY = element.coords.node[n][1];
            SIMD coordsZ = element.coords.node[n][2];
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);
            SIMD dNZ = load1(element.dN[gp][n][2]);

            J[0] = J[0] + dNX * coordsX;
            J[3] = J[3] + dNX * coordsY;
            J[6] = J[6] + dNX * coordsZ;
            J[1] = J[1] + dNY * coordsX;
            J[4] = J[4] + dNY * coordsY;
            J[7] = J[7] + dNY * coordsZ;
            J[2] = J[2] + dNZ * coordsX;
            J[5] = J[5] + dNZ * coordsY;
            J[8] = J[8] + dNZ * coordsZ;
        }

        inv33(J, element.det, element.invJ);

        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);
            SIMD dNZ = load1(element.dN[gp][n][2]);
            element.dND[n * 3 + 0] = element.invJ[0] * dNX + element.invJ[3] * dNY + element.invJ[6] * dNZ;
            element.dND[n * 3 + 1] = element.invJ[1] * dNX + element.invJ[4] * dNY + element.invJ[7] * dNZ;
            element.dND[n * 3 + 2] = element.invJ[2] * dNX + element.invJ[5] * dNY + element.invJ[8] * dNZ;
        }
    }

    template <typename Element>
    void displacement(Element &element, size_t gp)
    {
        SIMD J[9];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD coordsZ = element.coords.node[n][2] + element.displacement[n][2];
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);
            SIMD dNZ = load1(element.dN[gp][n][2]);

            J[0] = J[0] + dNX * coordsX;
            J[3] = J[3] + dNX * coordsY;
            J[6] = J[6] + dNX * coordsZ;
            J[1] = J[1] + dNY * coordsX;
            J[4] = J[4] + dNY * coordsY;
            J[7] = J[7] + dNY * coordsZ;
            J[2] = J[2] + dNZ * coordsX;
            J[5] = J[5] + dNZ * coordsY;
            J[8] = J[8] + dNZ * coordsZ;
        }

        inv33(J, element.det, element.invJ);

        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);
            SIMD dNZ = load1(element.dN[gp][n][2]);
            element.dND[n * 3 + 0] = element.invJ[0] * dNX + element.invJ[3] * dNY + element.invJ[6] * dNZ;
            element.dND[n * 3 + 1] = element.invJ[1] * dNX + element.invJ[4] * dNY + element.invJ[7] * dNZ;
            element.dND[n * 3 + 2] = element.invJ[2] * dNX + element.invJ[5] * dNY + element.invJ[8] * dNZ;
        }
    }

    template <typename Element>
    void displacementInNodes(Element &element, size_t node)
    {
        SIMD J[9];
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD coordsZ = element.coords.node[n][2] + element.displacement[n][2];
            SIMD dNX = load1(element.dNN[node][n][0]);
            SIMD dNY = load1(element.dNN[node][n][1]);
            SIMD dNZ = load1(element.dNN[node][n][2]);

            J[0] = J[0] + dNX * coordsX;
            J[3] = J[3] + dNX * coordsY;
            J[6] = J[6] + dNX * coordsZ;
            J[1] = J[1] + dNY * coordsX;
            J[4] = J[4] + dNY * coordsY;
            J[7] = J[7] + dNY * coordsZ;
            J[2] = J[2] + dNZ * coordsX;
            J[5] = J[5] + dNZ * coordsY;
            J[8] = J[8] + dNZ * coordsZ;
        }

        inv33(J, element.det, element.invJ);

        for (size_t n = 0; n < nodes; ++n) {
            SIMD dNX = load1(element.dNN[node][n][0]);
            SIMD dNY = load1(element.dNN[node][n][1]);
            SIMD dNZ = load1(element.dNN[node][n][2]);
            element.dND[n * 3 + 0] = element.invJ[0] * dNX + element.invJ[3] * dNY + element.invJ[6] * dNZ;
            element.dND[n * 3 + 1] = element.invJ[1] * dNX + element.invJ[4] * dNY + element.invJ[7] * dNZ;
            element.dND[n * 3 + 2] = element.invJ[2] * dNX + element.invJ[5] * dNY + element.invJ[8] * dNZ;
        }
    }
};

template <size_t nodes>
struct IntegrationKernel<nodes, 3, 2>: Integration {
    IntegrationKernel(const Integration &base): Integration(base) {}

    template <typename Element>
    void coords(Element &element, size_t gp)
    {
        SIMD dND0, dND1, dND2, dND3, dND4, dND5;
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0];
            SIMD coordsY = element.coords.node[n][1];
            SIMD coordsZ = element.coords.node[n][2];
            SIMD dNX = load1(element.dN[gp][n][0]);
            SIMD dNY = load1(element.dN[gp][n][1]);

            dND0 = dND0 + dNX * coordsX;
            dND1 = dND1 + dNX * coordsY;
            dND2 = dND2 + dNX * coordsZ;
            dND3 = dND3 + dNY * coordsX;
            dND4 = dND4 + dNY * coordsY;
            dND5 = dND5 + dNY * coordsZ;
        }

        SIMD x = dND1 * dND5 - dND2 * dND4;
        SIMD y = dND2 * dND3 - dND0 * dND5;
        SIMD z = dND0 * dND4 - dND1 * dND3;
        SIMD res = x * x + y * y + z * z;
        for (size_t s = 0; s < SIMD::size; ++s) {
            element.det[s] = std::sqrt(res[s]);
        }
        element.normal[0] = x / element.det;
        element.normal[1] = y / element.det;
        element.normal[2] = z / element.det;
    }
};

template <size_t nodes>
struct IntegrationKernel<nodes, 3, 1>: Integration {
    IntegrationKernel(const Integration &base): Integration(base) {}

    template <typename Element>
    void coords(Element &element, size_t gp)
    {
        SIMD dND0, dND1, dND2;
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0];
            SIMD coordsY = element.coords.node[n][1];
            SIMD coordsZ = element.coords.node[n][2];
            SIMD dNX = load1(element.dN[gp][n][0]);

            dND0 = dND0 + dNX * coordsX;
            dND1 = dND1 + dNX * coordsY;
            dND2 = dND2 + dNX * coordsZ;
        }

        SIMD res = dND0 * dND0 + dND1 * dND1 + dND2 * dND2;
        for (size_t s = 0; s < SIMD::size; ++s) {
            element.det[s] = std::sqrt(res[s]);
        }
    }
};

template <size_t nodes>
struct IntegrationKernel<nodes, 2, 1>: Integration {
    IntegrationKernel(const Integration &base): Integration(base) {}

    template <typename Element>
    void coords(Element &element, size_t gp)
    {
        SIMD dND0, dND1;
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0];
            SIMD coordsY = element.coords.node[n][1];
            SIMD dNX = load1(element.dN[gp][n][0]);

            dND0 = dND0 + dNX * coordsX;
            dND1 = dND1 + dNX * coordsY;
        }

        SIMD res = dND0 * dND0 + dND1 * dND1;
        for (size_t s = 0; s < SIMD::size; ++s) {
            element.det[s] = std::sqrt(res[s]);
        }
        element.normal[0] = -dND1 / element.det;
        element.normal[1] =  dND0 / element.det;
    }
};

}

#endif /* SRC_PHYSICS_ASSEMBLER_OPERATORS_INTEGRATION_H_ */
