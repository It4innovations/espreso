
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_NORMAL_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_NORMAL_H_

#include "element.h"
#include "subkernel.h"
#include "basis/containers/serializededata.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"

namespace espreso {

struct Normal: SubKernel {
    const char* name() const { return "Normal"; }

    serializededata<esint, esint>::const_iterator enodes, end;
    double *target, *multiplicity;

    Normal()
    : enodes(info::mesh->elements->nodes->cbegin()),
      end(info::mesh->elements->nodes->cend()),
      target(nullptr),
      multiplicity(nullptr)
    {
        isconst = false;
        action = SubKernel::PREPROCESS | SubKernel::ASSEMBLE | SubKernel::REASSEMBLE;
    }

    void activate(serializededata<esint, esint>::const_iterator enodes, serializededata<esint, esint>::const_iterator end, double *target, double *multiplicity)
    {
        this->enodes = enodes;
        this->end = end;
        this->target = target;
        this->multiplicity = multiplicity;
        this->isactive = 1;
    }
};

template <size_t nodes, size_t ndim, size_t edim> struct NormalKernel;

template <size_t nodes>
struct NormalKernel<nodes, 2, 1>: Normal {
    NormalKernel(const Normal &base): Normal(base) { }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD dND0, dND1;
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD dNX = load1(element.dN[gp][n][0]);

            dND0 = dND0 + dNX * coordsX;
            dND1 = dND1 + dNX * coordsY;
        }

        SIMD res = rsqrt14(dND0 * dND0 + dND1 * dND1);
        element.normal[0] = -dND1 * res;
        element.normal[1] =  dND0 * res;
    }

    template <typename Element>
    void simd(Element &element)
    {
        // TODO: compute normal in nodes
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t n = 0; n < nodes; ++n) {
                for (size_t d = 0; d < 2; ++d) {
                    target[3 * enodes->at(n) + d] += element.normal[d][s] * multiplicity[enodes->at(n)];
                }
            }
        }
    }
};

template <size_t nodes>
struct NormalKernel<nodes, 3, 2>: Normal {
    NormalKernel(const Normal &base): Normal(base) { }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {
        SIMD dND0, dND1, dND2, dND3, dND4, dND5;
        for (size_t n = 0; n < nodes; ++n) {
            SIMD coordsX = element.coords.node[n][0] + element.displacement[n][0];
            SIMD coordsY = element.coords.node[n][1] + element.displacement[n][1];
            SIMD coordsZ = element.coords.node[n][2] + element.displacement[n][2];
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
        SIMD res = rsqrt14(x * x + y * y + z * z);
        element.normal[0] = x * res;
        element.normal[1] = y * res;
        element.normal[2] = z * res;
    }

    template <typename Element>
    void simd(Element &element)
    {
        // TODO: compute normal in nodes
        for (size_t s = 0; s < SIMD::size; ++s, ++enodes) {
            if (enodes == end) break;
            for (size_t n = 0; n < nodes; ++n) {
                for (size_t d = 0; d < 3; ++d) {
                    target[3 * enodes->at(n) + d] += element.normal[d][s] * multiplicity[enodes->at(n)];
                }
            }
        }
    }
};

template <size_t nodes>
struct NormalKernel<nodes, 3, 1>: Normal {
    NormalKernel(const Normal &base): Normal(base) { }

    template <typename Element>
    void simd(Element &element, size_t gp)
    {

    }

    template <typename Element>
    void simd(Element &element)
    {

    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_OP_NORMAL_H_ */
