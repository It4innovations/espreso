
#ifndef SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_ELEMENT_H_
#define SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_ELEMENT_H_

#include "analysis/assembler/general/element.h"
#include "analysis/assembler/general/subkernel.h"
#include "analysis/assembler/structuralmechanics/operators.h"
#include "esinfo/stepinfo.h"
#include "mesh/element.h"

namespace espreso {

struct StructuralMechanicsGPC {
    enum: int {
        POINT1    =  1,

        LINE2     =  2,
        TRIANGLE3 =  4,
        SQUARE4   =  4,
        TETRA4    =  8,
        PYRAMID5  =  8,
        PRISMA6   =  8,
        HEXA8     =  8,

        LINE3     =  4,
        TRIANGLE6 =  6,
        SQUARE8   = 16,
        TETRA10   = 15,
        PYRAMID13 = 27,
        PRISMA15  =  9,
        HEXA20    = 27,
    };
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct StructuralMechanicsElement: public GeneralElement<nodes, gps, ndim, edim> {
    struct {
        alignas(SIMD::size * sizeof(double)) SIMD density;
        alignas(SIMD::size * sizeof(double)) SIMD heatCapacity;

        // elasticity
        alignas(SIMD::size * sizeof(double)) SIMD youngModulus[ndim];
        alignas(SIMD::size * sizeof(double)) SIMD poissonRatio[ndim];
        alignas(SIMD::size * sizeof(double)) SIMD shearModulus[ndim];

        // plasticity
        alignas(SIMD::size * sizeof(double)) SIMD initialYieldStress;
        alignas(SIMD::size * sizeof(double)) SIMD isotropicHardening;
        alignas(SIMD::size * sizeof(double)) SIMD kinematicHardening;
        alignas(SIMD::size * sizeof(double)) SIMD sigma;

        alignas(SIMD::size * sizeof(double)) SIMD acceleration   [ndim];
        alignas(SIMD::size * sizeof(double)) SIMD angularVelocity[3];
    } ecf;


    alignas(SIMD::size * sizeof(double)) SIMD F[ndim * ndim];
    alignas(SIMD::size * sizeof(double)) SIMD C2[ndim * ndim]; // 2D = 16, 3D = 36
    // PLANE = 9, AXISYMMETRIC = 16, VOLUME = 36
    alignas(SIMD::size * sizeof(double)) SIMD vC4[ndim * ndim * 4]; // 2D = 16, 3D = 36
    alignas(SIMD::size * sizeof(double)) SIMD vS[ndim + ndim * (ndim - 1) / 2];

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD center[ndim];
        alignas(SIMD::size * sizeof(double)) SIMD cossin[ndim * 4];
    } rotation;

    alignas(SIMD::size * sizeof(double)) SIMD K[ndim * nodes * ndim * nodes];
    alignas(SIMD::size * sizeof(double)) SIMD M[ndim * nodes * ndim * nodes];
    alignas(SIMD::size * sizeof(double)) SIMD C[ndim * nodes * ndim * nodes];
    alignas(SIMD::size * sizeof(double)) SIMD f[ndim * nodes], imf[ndim * nodes];
    alignas(SIMD::size * sizeof(double)) SIMD nf[ndim * nodes];
    alignas(SIMD::size * sizeof(double)) SIMD stress[nodes * 13];

    alignas(SIMD::size * sizeof(double)) SIMD displacement[nodes][ndim];
    alignas(SIMD::size * sizeof(double)) SIMD velocity[nodes][ndim];

    StructuralMechanicsElement()
    {
        for (size_t i = 0; i < ndim * nodes * ndim * nodes; ++i) {
            K[i] = zeros();
            M[i] = zeros();
            C[i] = zeros();
        }
    }
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct StructuralMechanicsBoundary: public GeneralBoundary<nodes, gps, ndim, edim> {
    struct {
        alignas(SIMD::size * sizeof(double)) SIMD normalPressure;
        struct {
            alignas(SIMD::size * sizeof(double)) SIMD pressure;
            alignas(SIMD::size * sizeof(double)) SIMD direction[ndim];
        } pressure;
    } ecf;

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD pressure[nodes];
        alignas(SIMD::size * sizeof(double)) SIMD stress[nodes][ndim];
    } coupling;

    alignas(SIMD::size * sizeof(double)) SIMD displacement[nodes][ndim];
    alignas(SIMD::size * sizeof(double)) SIMD normal[ndim];
    alignas(SIMD::size * sizeof(double)) SIMD f[ndim * nodes];
};

template <size_t ndim> struct StructuralMechanicsDirichlet: public GeneralNode<ndim> {
    struct {
        alignas(SIMD::size * sizeof(double)) SIMD force[ndim];
        alignas(SIMD::size * sizeof(double)) SIMD harmonicForceMag[ndim];
        alignas(SIMD::size * sizeof(double)) SIMD harmonicForceCos[ndim];
        alignas(SIMD::size * sizeof(double)) SIMD harmonicForceSin[ndim];
    } ecf;

    struct {
        alignas(SIMD::size * sizeof(double)) SIMD node[ndim];
    } displacement;

    alignas(SIMD::size * sizeof(double)) SIMD f[ndim], imf[ndim];
};

template <Element::CODE code> void runElement(const step::Step &step, StructuralMechanicsElementOperators &operators, SubKernel::Action action);
template <Element::CODE code> void runBoundary(const step::Step &step, StructuralMechanicsFaceOperators &operators, SubKernel::Action action);
template <Element::CODE code> void runNode(const step::Step &step, StructuralMechanicsNodeOperators &operators, SubKernel::Action action);

}

#endif /* SRC_ANALYSIS_ASSEMBLER_STRUCTURALMECHANICS_ELEMENT_H_ */

