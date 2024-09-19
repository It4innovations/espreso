
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BASIS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BASIS_H_

#include "element.h"
#include "subkernel.h"
#include "mesh/element.h"

namespace espreso {

struct Basis: SubKernel {
    const char* name() const { return "Basis"; }

    Basis()
    {
        action = SubKernel::ASSEMBLE | SubKernel::REASSEMBLE | SubKernel::ITERATION | SubKernel::SOLUTION;
    }
};

template <Element::CODE code, size_t nodes, size_t gps, size_t edim> struct GaussPoints;

template <Element::CODE code, size_t nodes, size_t gps, size_t edim>
struct BasisKernel: Basis {
    BasisKernel(const Basis &base): Basis(base) { }

    template <typename Element>
    void simd(Element &element)
    {
        for (size_t gp = 0; gp < gps; ++gp) {
            element.w[gp] = GaussPoints<code, nodes, gps, edim>::w[gp];
            for (size_t n = 0; n < nodes; ++n) {
                element.N[gp][n] = GaussPoints<code, nodes, gps, edim>::N[gp * nodes + n];
                for (size_t d = 0; d < edim; ++d) {
                    element.dN[gp][n][d] = GaussPoints<code, nodes, gps, edim>::dN[gp * edim * nodes + d * nodes + n];
                }
            }
        }
    }

    static void setCenter(double &w, double N[], double dN[][edim])
    {
        w = GaussPoints<code, nodes, gps, edim>::cw;
        for (size_t n = 0; n < nodes; ++n) {
            N[n] = GaussPoints<code, nodes, gps, edim>::cN[n];
            for (size_t d = 0; d < edim; ++d) {
                dN[n][d] = GaussPoints<code, nodes, gps, edim>::cdN[d * nodes + n];
            }
        }
    }
};



template <Element::CODE code, size_t nodes, size_t gps, size_t edim>
struct BasisCenterKernel {
    static void set(double *w, double *N, double **dN)
    {

    }
};

template<>
struct GaussPoints<Element::CODE::LINE2, 2, 2, 1> {

    constexpr static int nodes = 2, gps = 2, edim = 1;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double s)
    {
        _N[0] = (1 - s) * 0.5;
        _N[1] = (1 + s) * 0.5;

        _dN[0] = -0.5;
        _dN[1] =  0.5;
    }

    static void set()
    {
        double s[2] = { 1 / sqrt(3), -1 / sqrt(3) };
        for (int gp = 0; gp < gps; gp++) {
            w[gp] = 1;
            set(N + gp * nodes, dN + edim * gp * nodes, s[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2> {

    constexpr static int nodes = 3, gps = 6, edim = 2;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double s, double t)
    {
        _N[0] = 1 - s - t;
        _N[1] = s;
        _N[2] = t;

        _dN[0 * nodes + 0] = -1;
        _dN[0 * nodes + 1] =  1;
        _dN[0 * nodes + 2] =  0;

        _dN[1 * nodes + 0] = -1;
        _dN[1 * nodes + 1] =  0;
        _dN[1 * nodes + 2] =  1;
    }

    static void set()
    {
        double s[6] = { 0.445948490915965, 0.445948490915965, 0.108103018168070, 0.091576213509771, 0.091576213509771, 0.816847572980459 };
        double t[6] = { 0.445948490915965, 0.108103018168070, 0.445948490915965, 0.091576213509771, 0.816847572980459, 0.091576213509771 };

        w[0] = 0.111690794839005;
        w[1] = 0.111690794839005;
        w[2] = 0.111690794839005;
        w[3] = 0.054975871827661;
        w[4] = 0.054975871827661;
        w[5] = 0.054975871827661;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::SQUARE4, 4, 4, 2> {

    constexpr static int nodes = 4, gps = 4, edim = 2;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double s, double t)
    {
        _N[0] = 0.25 * (1 - s) * (1 - t);
        _N[1] = 0.25 * (s + 1) * (1 - t);
        _N[2] = 0.25 * (s + 1) * (t + 1);
        _N[3] = 0.25 * (1 - s) * (t + 1);

        _dN[0 * nodes + 0] = 0.25 * ( t - 1);
        _dN[0 * nodes + 1] = 0.25 * (-t + 1);
        _dN[0 * nodes + 2] = 0.25 * ( t + 1);
        _dN[0 * nodes + 3] = 0.25 * (-t - 1);

        _dN[1 * nodes + 0] = 0.25 * ( s - 1);
        _dN[1 * nodes + 1] = 0.25 * (-s - 1);
        _dN[1 * nodes + 2] = 0.25 * ( s + 1);
        _dN[1 * nodes + 3] = 0.25 * (-s + 1);
    }

    static void set()
    {
        double CsQ_scale = 0.577350269189626;
        double s[4] = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
        double t[4] = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };

        for (int gp = 0; gp < gps; gp++) {
            w[gp] = 1;
            set(N + gp * nodes, dN + edim * gp * nodes, s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::TETRA4, 4, 4, 3> {

    constexpr static int nodes = 4, gps = 4, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[0] = r;
        _N[1] = s;
        _N[2] = t;
        _N[3] = 1.0 - r - s - t;

        _dN[0 * nodes + 0] =  1.0;
        _dN[0 * nodes + 1] =  0.0;
        _dN[0 * nodes + 2] =  0.0;
        _dN[0 * nodes + 3] = -1.0;

        _dN[1 * nodes + 0] =  0.0;
        _dN[1 * nodes + 1] =  0.0;
        _dN[1 * nodes + 2] =  1.0;
        _dN[1 * nodes + 3] = -1.0;

        _dN[2 * nodes + 0] =  0.0;
        _dN[2 * nodes + 1] =  1.0;
        _dN[2 * nodes + 2] =  0.0;
        _dN[2 * nodes + 3] = -1.0;
    }

    static void set()
    {
        double r[4] = { 0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105 };
        double s[4] = { 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685 };
        double t[4] = { 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105 };

        for (int gp = 0; gp < gps; gp++) {
            w[gp] = 1.0 / 24.0;
            set(N + gp * nodes, dN + edim * gp * nodes, r[gp], s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3> {

    constexpr static int nodes = 5, gps = 8, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[0] = 0.125 * ((1 - r) * (1 - s) * (1 - t));
        _N[1] = 0.125 * ((1 + r) * (1 - s) * (1 - t));
        _N[2] = 0.125 * ((1 + r) * (1 + s) * (1 - t));
        _N[3] = 0.125 * ((1 - r) * (1 + s) * (1 - t));
        _N[4] = 0.125 * ( 4 * (1 + t));

        _dN[0 * nodes + 0] = 0.125 * (-(1. - s) * (1. - t));
        _dN[0 * nodes + 1] = 0.125 * ( (1. - s) * (1. - t));
        _dN[0 * nodes + 2] = 0.125 * ( (1. + s) * (1. - t));
        _dN[0 * nodes + 3] = 0.125 * (-(1. + s) * (1. - t));
        _dN[0 * nodes + 4] = 0;

        _dN[1 * nodes + 0] = 0.125 * (-(1. - r) * (1. - t));
        _dN[1 * nodes + 1] = 0.125 * (-(1. + r) * (1. - t));
        _dN[1 * nodes + 2] = 0.125 * ( (1. + r) * (1. - t));
        _dN[1 * nodes + 3] = 0.125 * ( (1. - r) * (1. - t));
        _dN[1 * nodes + 4] = 0;

        _dN[2 * nodes + 0] = 0.125 * (-(1. - r) * (1. - s));
        _dN[2 * nodes + 1] = 0.125 * (-(1. + r) * (1. - s));
        _dN[2 * nodes + 2] = 0.125 * (-(1. + r) * (1. + s));
        _dN[2 * nodes + 3] = 0.125 * (-(1. - r) * (1. + s));
        _dN[2 * nodes + 4] = 0.125 * (4.0);
    }

    static void set()
    {
        double v = 0.577350269189625953;
        double r[8] = {  v,  v,  v,  v, -v, -v, -v, -v };
        double s[8] = { -v, -v,  v,  v, -v, -v,  v,  v };
        double t[8] = { -v,  v, -v,  v, -v,  v, -v,  v };

        for (int gp = 0; gp < gps; gp++) {
            w[gp] = 1;
            set(N + gp * nodes, dN + edim * gp * nodes, r[gp], s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::PRISMA6, 6, 9, 3> {

    constexpr static int nodes = 6, gps = 9, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[0] = 0.5 * ((1.0 - t) * (1.0 - r - s));
        _N[1] = 0.5 * ((1.0 - t) * r);
        _N[2] = 0.5 * ((1.0 - t) * s);
        _N[3] = 0.5 * ((1.0 + t) * (1.0 - r - s));
        _N[4] = 0.5 * ((1.0 + t) * r);
        _N[5] = 0.5 * ((1.0 + t) * s);

        _dN[0 * nodes + 0] =  t / 2.0 - 1.0 / 2.0;
        _dN[0 * nodes + 1] = -t / 2.0 + 1.0 / 2.0;
        _dN[0 * nodes + 2] =  0.0;
        _dN[0 * nodes + 3] = -t / 2.0 - 1.0 / 2.0;
        _dN[0 * nodes + 4] =  t / 2.0 + 1.0 / 2.0;
        _dN[0 * nodes + 5] =  0;

        _dN[1 * nodes + 0] =  t / 2.0 - 1.0 / 2.0;
        _dN[1 * nodes + 1] =  0.0;
        _dN[1 * nodes + 2] = -t / 2.0 + 1.0 / 2.0;
        _dN[1 * nodes + 3] = -t / 2.0 - 1.0 / 2.0;
        _dN[1 * nodes + 4] =  0.0;
        _dN[1 * nodes + 5] =  t / 2.0 + 1.0 / 2.0;

        _dN[2 * nodes + 0] =  r / 2.0 + s / 2.0 - 1.0 / 2.0;
        _dN[2 * nodes + 1] = -r / 2.0;
        _dN[2 * nodes + 2] =          - s / 2.0;
        _dN[2 * nodes + 3] = -r / 2.0 - s / 2.0 + 1.0 / 2.0;
        _dN[2 * nodes + 4] =  r / 2.0;
        _dN[2 * nodes + 5] =            s / 2.0;
    }

    static void set()
    {
        double v1 = 1.0 / 6.0;
        double v2 = 4.0 / 6.0;
        double v3 = sqrt(3.0 / 5.0);
        double v4 = 0.0;
        double r[9] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
        double s[9] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
        double t[9] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };

        w[0] = 5.0 / 54.0;
        w[1] = 5.0 / 54.0;
        w[2] = 5.0 / 54.0;
        w[3] = 8.0 / 54.0;
        w[4] = 8.0 / 54.0;
        w[5] = 8.0 / 54.0;
        w[6] = 5.0 / 54.0;
        w[7] = 5.0 / 54.0;
        w[8] = 5.0 / 54.0;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, r[gp], s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::HEXA8, 8, 8, 3> {

    constexpr static int nodes = 8, gps = 8, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[0] = 0.125 * (1 - r) * (1 - s) * (1 - t);
        _N[1] = 0.125 * (r + 1) * (1 - s) * (1 - t);
        _N[2] = 0.125 * (r + 1) * (s + 1) * (1 - t);
        _N[3] = 0.125 * (1 - r) * (s + 1) * (1 - t);
        _N[4] = 0.125 * (1 - r) * (1 - s) * (t + 1);
        _N[5] = 0.125 * (r + 1) * (1 - s) * (t + 1);
        _N[6] = 0.125 * (r + 1) * (s + 1) * (t + 1);
        _N[7] = 0.125 * (1 - r) * (s + 1) * (t + 1);

        _dN[0 * nodes + 0] = 0.125 * (-(1 - s) * (1 - t));
        _dN[0 * nodes + 1] = 0.125 * ( (1 - s) * (1 - t));
        _dN[0 * nodes + 2] = 0.125 * ( (1 + s) * (1 - t));
        _dN[0 * nodes + 3] = 0.125 * (-(1 + s) * (1 - t));
        _dN[0 * nodes + 4] = 0.125 * (-(1 - s) * (1 + t));
        _dN[0 * nodes + 5] = 0.125 * ( (1 - s) * (1 + t));
        _dN[0 * nodes + 6] = 0.125 * ( (1 + s) * (1 + t));
        _dN[0 * nodes + 7] = 0.125 * (-(1 + s) * (1 + t));

        _dN[1 * nodes + 0] = 0.125 * (-(1 - r) * (1 - t));
        _dN[1 * nodes + 1] = 0.125 * (-(1 + r) * (1 - t));
        _dN[1 * nodes + 2] = 0.125 * ( (1 + r) * (1 - t));
        _dN[1 * nodes + 3] = 0.125 * ( (1 - r) * (1 - t));
        _dN[1 * nodes + 4] = 0.125 * (-(1 - r) * (1 + t));
        _dN[1 * nodes + 5] = 0.125 * (-(1 + r) * (1 + t));
        _dN[1 * nodes + 6] = 0.125 * ( (1 + r) * (1 + t));
        _dN[1 * nodes + 7] = 0.125 * ( (1 - r) * (1 + t));

        _dN[2 * nodes + 0] = 0.125 * (-(1 - r) * (1 - s));
        _dN[2 * nodes + 1] = 0.125 * (-(1 + r) * (1 - s));
        _dN[2 * nodes + 2] = 0.125 * (-(1 + r) * (1 + s));
        _dN[2 * nodes + 3] = 0.125 * (-(1 - r) * (1 + s));
        _dN[2 * nodes + 4] = 0.125 * ( (1 - r) * (1 - s));
        _dN[2 * nodes + 5] = 0.125 * ( (1 + r) * (1 - s));
        _dN[2 * nodes + 6] = 0.125 * ( (1 + r) * (1 + s));
        _dN[2 * nodes + 7] = 0.125 * ( (1 - r) * (1 + s));
    }

    static void set()
    {
        double CsQ_scale = 1 / std::sqrt(3);

        for (int gp = 0; gp < gps; gp++) {
            double r = (gp & 4) ? CsQ_scale : -CsQ_scale;
            double s = (gp & 2) ? CsQ_scale : -CsQ_scale;
            double t = (gp & 1) ? CsQ_scale : -CsQ_scale;

            w[gp] = 1;
            set(N + gp * nodes, dN + edim * gp * nodes, r, s, t);
        }

        { // center
            cw = 1;
            set(cN, cdN, 0, 0, 0);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::LINE3, 3, 3, 1> {

    constexpr static int nodes = 3, gps = 3, edim = 1;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double s)
    {
        _N[0] = 0.5 * (s - 1) * s;
        _N[1] = 0.5 * (s + 1) * s;
        _N[2] = 1 - s * s;

        _dN[0] = s - 0.5;
        _dN[1] = s + 0.5;
        _dN[2] = -2 * s;
    }

    static void set()
    {
        double s[3] = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

        w[0] = 5/9.0;
        w[1] = 8/9.0;
        w[2] = 5/9.0;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, s[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2> {

    constexpr static int nodes = 6, gps = 6, edim = 2;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double s, double t)
    {
        _N[0] = (1.0 - s - t) * (1.0 - 2.0 * (s + t));
        _N[1] = -(s) * (1.0 - 2.0 * s);
        _N[2] = -(t) * (1.0 - 2.0 * t);
        _N[3] = 4.0 * (s) * (1.0 - s - t);
        _N[4] = 4.0 * (s) * (t);
        _N[5] = 4.0 * (t) * (1.0 - s - t);

        _dN[0 * nodes + 0] = -3.0 + 4.0 * s + 4.0 * t;
        _dN[0 * nodes + 1] = -1.0 + 4.0 * s;
        _dN[0 * nodes + 2] = 0.0;
        _dN[0 * nodes + 3] = 4.0 - 8.0 * s - 4.0 * t;
        _dN[0 * nodes + 4] = 4.0 * t;
        _dN[0 * nodes + 5] = -4.0 * t;

        _dN[1 * nodes + 0] = -3.0 + 4.0 * s + 4.0 * t;
        _dN[1 * nodes + 1] = 0.0;
        _dN[1 * nodes + 2] = -1.0 + 4.0 * t;
        _dN[1 * nodes + 3] = -4.0 * s;
        _dN[1 * nodes + 4] = 4.0 * s;
        _dN[1 * nodes + 5] = 4.0 - 4.0 * s - 8.0 * t;
    }

    static void set()
    {
        double s[6] = { 0.091576213509771, 0.816847572980459, 0.091576213509771, 0.445948490915965, 0.108103018168070, 0.445948490915965 };
        double t[6] = { 0.091576213509771, 0.091576213509771, 0.816847572980459, 0.445948490915965, 0.445948490915965, 0.108103018168070 };

        w[0] = 0.109951743655322 / 2.0;
        w[1] = 0.109951743655322 / 2.0;
        w[2] = 0.109951743655322 / 2.0;
        w[3] = 0.223381589678011 / 2.0;
        w[4] = 0.223381589678011 / 2.0;
        w[5] = 0.223381589678011 / 2.0;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::SQUARE8, 8, 9, 2> {

    constexpr static int nodes = 8, gps = 9, edim = 2;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double s, double t)
    {
        _N[0] = -.25 * (s - 1) * (t - 1) * (s + t + 1);
        _N[1] =  .25 * (t - 1) * (-s * s + t * s + t + 1);
        _N[2] =  .25 * (s + 1) * (t + 1) * (s + t - 1);
        _N[3] =  .25 * (s - 1) * (s - t + 1) * (t + 1);
        _N[4] =  .5  * (s * s - 1) * (t - 1);
        _N[5] = -.5  * (s + 1) * (t * t - 1);
        _N[6] = -.5  * (s * s - 1) * (t + 1);
        _N[7] =  .5  * (s - 1) * (t * t - 1);

        _dN[0 * nodes + 0] = -((2 * s + t) * (t - 1)) * .25;
        _dN[0 * nodes + 1] = -((2 * s - t) * (t - 1)) * .25;
        _dN[0 * nodes + 2] =  ((2 * s + t) * (t + 1)) * .25;
        _dN[0 * nodes + 3] =  ((2 * s - t) * (t + 1)) * .25;
        _dN[0 * nodes + 4] = s * (t - 1);
        _dN[0 * nodes + 5] = .5 - t * t * .5;
        _dN[0 * nodes + 6] = -s * (t + 1);
        _dN[0 * nodes + 7] = t * t * .5 - .5;

        _dN[1 * nodes + 0] = -((s + 2 * t) * (s - 1)) * .25;
        _dN[1 * nodes + 1] = -((s - 2 * t) * (s + 1)) * .25;
        _dN[1 * nodes + 2] =  ((s + 2 * t) * (s + 1)) * .25;
        _dN[1 * nodes + 3] =  ((s - 2 * t) * (s - 1)) * .25;
        _dN[1 * nodes + 4] = s * s * .5 - .5;
        _dN[1 * nodes + 5] = -t * (s + 1);
        _dN[1 * nodes + 6] = .5 - s * s * .5;
        _dN[1 * nodes + 7] = t * (s - 1);
    }

    static void set()
    {
        double v = sqrt(0.6);
        double s[9] = { -v,  v,  v, -v,  0,  v,  0, -v, 0 };
        double t[9] = { -v, -v,  v,  v, -v,  0,  v,  0, 0 };

        w[0] = 25.0 / 81.0;
        w[1] = 25.0 / 81.0;
        w[2] = 25.0 / 81.0;
        w[3] = 25.0 / 81.0;
        w[4] = 40.0 / 81.0;
        w[5] = 40.0 / 81.0;
        w[6] = 40.0 / 81.0;
        w[7] = 40.0 / 81.0;
        w[8] = 64.0 / 81.0;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::TETRA10, 10, 15, 3> {

    constexpr static int nodes = 10, gps = 15, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[0] = r * (2.0 * r - 1.0);
        _N[1] = s * (2.0 * s - 1.0);
        _N[2] = t * (2.0 * t - 1.0);
        _N[3] = 2.0 * r * r + 4.0 * r * s + 4.0 * r * t - 3.0 * r + 2.0* s * s + 4.0 * s * t - 3.0 * s + 2.0 * t * t - 3.0 * t + 1.0;
        _N[4] = 4.0 * r * s;
        _N[5] = 4.0 * s * t;
        _N[6] = 4.0 * r * t;
        _N[7] = r * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
        _N[8] = s * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
        _N[9] = t * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);

        _dN[0 * nodes + 0] = 4.0 * r - 1.0;
        _dN[0 * nodes + 1] = 0;
        _dN[0 * nodes + 2] = 0;
        _dN[0 * nodes + 3] = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
        _dN[0 * nodes + 4] = 4.0 * s;
        _dN[0 * nodes + 5] = 0;
        _dN[0 * nodes + 6] = 4.0 * t;
        _dN[0 * nodes + 7] = -8.0 * r - 4.0 * s - 4.0 * t + 4.0;
        _dN[0 * nodes + 8] = -4.0 * s;
        _dN[0 * nodes + 9] = -4.0 * t;

        _dN[1 * nodes + 0] = 0;
        _dN[1 * nodes + 1] = 0;
        _dN[1 * nodes + 2] = 4.0 * t - 1.0;
        _dN[1 * nodes + 3] = 4.0 * r + 4.0 * s + 4.0* t  - 3.0;
        _dN[1 * nodes + 4] = 0;
        _dN[1 * nodes + 5] = 4.0 * s;
        _dN[1 * nodes + 6] = 4.0 * r;
        _dN[1 * nodes + 7] = -4.0 * r;
        _dN[1 * nodes + 8] = -4.0 * s;
        _dN[1 * nodes + 9] = -4.0 * r - 4.0 * s - 8.0 * t + 4.0;

        _dN[2 * nodes + 0] = 0;
        _dN[2 * nodes + 1] = 4.0 * s - 1.0;
        _dN[2 * nodes + 2] = 0 ;
        _dN[2 * nodes + 3] = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
        _dN[2 * nodes + 4] = 4.0 * r;
        _dN[2 * nodes + 5] = 4.0 * t;
        _dN[2 * nodes + 6] = 0;
        _dN[2 * nodes + 7] = -4.0 * r;
        _dN[2 * nodes + 8] = -4.0 * r - 8.0 * s - 4.0 * t + 4.0;
        _dN[2 * nodes + 9] = -4.0 * t;
    }

    static void set()
    {
        double r[15] = {
                0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
                0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
                0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
                0.0665501535736643, 0.4334498464263357, 0.4334498464263357 };
        double s[15] = {
                0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
                0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
                0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
                0.4334498464263357, 0.0665501535736643, 0.4334498464263357 };
        double t[15] = {
                0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
                0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
                0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
                0.4334498464263357, 0.4334498464263357, 0.0665501535736643 };

        w[ 0] = 0.030283678097089;
        w[ 1] = 0.006026785714286;
        w[ 2] = 0.006026785714286;
        w[ 3] = 0.006026785714286;
        w[ 4] = 0.006026785714286;
        w[ 5] = 0.011645249086029;
        w[ 6] = 0.011645249086029;
        w[ 7] = 0.011645249086029;
        w[ 8] = 0.011645249086029;
        w[ 9] = 0.010949141561386;
        w[10] = 0.010949141561386;
        w[11] = 0.010949141561386;
        w[12] = 0.010949141561386;
        w[13] = 0.010949141561386;
        w[14] = 0.010949141561386;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, r[gp], s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3> {

    constexpr static int nodes = 13, gps = 14, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[ 0] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 - s) * (-1.0 - (0.5 * (1.0 - t)) * r - (0.5 * (1.0 - t)) * s));
        _N[ 1] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 - s) * (-1.0 + (0.5 * (1 - t))   * r - (0.5 * (1.0 - t)) * s));
        _N[ 2] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 + s) * (-1.0 + (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
        _N[ 3] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 + s) * (-1.0 - (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
        _N[ 4] = (1.0 - (0.5 * (1.0 - t))) * (1.0 - 2.0 * (0.5 * (1.0 - t)));
        _N[ 5] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - s) * (1.0 - r * r);
        _N[ 6] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + r) * (1.0 - s * s);
        _N[ 7] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + s) * (1.0 - r * r);
        _N[ 8] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - r) * (1.0 - s * s);
        _N[ 9] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r - s + r * s);
        _N[10] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r - s - r * s);
        _N[11] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r + s + r * s);
        _N[12] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r + s - r * s);

        _dN[0 * nodes +  0] = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
        _dN[0 * nodes +  1] = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
        _dN[0 * nodes +  2] =  (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0) + (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
        _dN[0 * nodes +  3] =  (t / 2.0 - 0.5) * (s + 1.0) * (r - 1.0) * (t / 8.0 - 1.0 / 8.0) - (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0);
        _dN[0 * nodes +  4] =  0.0;
        _dN[0 * nodes +  5] =  r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s - 1.0);
        _dN[0 * nodes +  6] = -((s * s - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
        _dN[0 * nodes +  7] = -r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s + 1.0);
        _dN[0 * nodes +  8] =  ((s * s - 1) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
        _dN[0 * nodes +  9] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
        _dN[0 * nodes + 10] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
        _dN[0 * nodes + 11] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);
        _dN[0 * nodes + 12] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);

        _dN[1 * nodes +  0] = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) - 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
        _dN[1 * nodes +  1] =  (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0) - (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) - s * (t / 2.0 - 1.0 / 2.0) + 1.0);
        _dN[1 * nodes +  2] =  (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) + 1.0) + (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
        _dN[1 * nodes +  3] = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s * (t / 2.0 - 1.0 / 2.0) - r * (t / 2.0 - 1.0 / 2.0) + 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0);
        _dN[1 * nodes +  4] =  0.0;
        _dN[1 * nodes +  5] =  ((r * r - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
        _dN[1 * nodes +  6] = -s * (t / 2.0 - 1.0 / 2.0) * (t / 2.0 - 0.5) * (r + 1.0);
        _dN[1 * nodes +  7] = -((r * r - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
        _dN[1 * nodes +  8] =  s * (t / 2.0 - 0.5) * (t / 2.0 - 0.5) * (r - 1.0);
        _dN[1 * nodes +  9] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);
        _dN[1 * nodes + 10] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
        _dN[1 * nodes + 11] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
        _dN[1 * nodes + 12] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);

        _dN[2 * nodes +  0] = -((r - 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0)) / 8.0 - (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
        _dN[2 * nodes +  1] = -((r + 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0)) / 8.0 - (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
        _dN[2 * nodes +  2] =  ((r + 1.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0)) / 8.0 + (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
        _dN[2 * nodes +  3] =  (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0) - ((r - 1.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0)) / 8.0;
        _dN[2 * nodes +  4] =  t + 0.5;
        _dN[2 * nodes +  5] =  ((r * r - 1.0) * (t / 2.0 - 0.5) * (s - 1.0)) / 2.0;
        _dN[2 * nodes +  6] = -((s * s - 1.0) * (t / 2.0 - 0.5) * (r + 1.0)) / 2.0;
        _dN[2 * nodes +  7] = -((r * r - 1.0) * (t / 2.0 - 0.5) * (s + 1.0)) / 2.0;
        _dN[2 * nodes +  8] =  ((s * s - 1.0) * (t / 2.0 - 0.5) * (r - 1.0)) / 2.0;
        _dN[2 * nodes +  9] =  ((t / 2.0 - 0.5) * (r + s - r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r + s - r * s - 1.0)) / 2.0;
        _dN[2 * nodes + 10] = -((t / 2.0 - 0.5) * (r - s - r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r - s - r * s + 1.0)) / 2.0;
        _dN[2 * nodes + 11] = -((t / 2.0 - 0.5) * (r + s + r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r + s + r * s + 1.0)) / 2.0;
        _dN[2 * nodes + 12] =  ((t / 2.0 - 0.5) * (r - s + r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r - s + r * s - 1.0)) / 2.0;
    }

    static void set()
    {
        double v1 = 0.758786910639329015;
        double v2 = 0.795822425754222018;
        double v3 = 0;
        double r[14] = { -v1,  v1,  v1, -v1, -v1,  v1,  v1, -v1,  v3,  v3,  v2, v3, -v2, v3 };
        double s[14] = { -v1, -v1,  v1,  v1, -v1, -v1,  v1,  v1,  v3, -v2,  v3, v2,  v3, v3 };
        double t[14] = { -v1, -v1, -v1, -v1,  v1,  v1,  v1,  v1, -v2,  v3,  v3, v3,  v3, v2 };

        w[ 0] = 0.335180055401662;
        w[ 1] = 0.335180055401662;
        w[ 2] = 0.335180055401662;
        w[ 3] = 0.335180055401662;
        w[ 4] = 0.335180055401662;
        w[ 5] = 0.335180055401662;
        w[ 6] = 0.335180055401662;
        w[ 7] = 0.335180055401662;
        w[ 8] = 0.886426592797784;
        w[ 9] = 0.886426592797784;
        w[10] = 0.886426592797784;
        w[11] = 0.886426592797784;
        w[12] = 0.886426592797784;
        w[13] = 0.886426592797784;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, r[gp], s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::PRISMA15, 15, 9, 3> {

    constexpr static int nodes = 15, gps = 9, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[ 0] = -(1.0 - r - s) * (1.0 - t) * (2.0 * r + 2.0 * s + t) / 2.0;
        _N[ 1] = r * (1.0 - t) * (2.0 * r - t - 2.0) / 2.0;
        _N[ 2] = s * (1.0 - t) * (2.0 * s - t - 2.0) / 2.0;
        _N[ 3] = -(1.0 - r - s) * (1.0 + t) * (2.0 * r + 2.0 * s - t) / 2.0;
        _N[ 4] = r * (t + 1.0) * (2.0 * r + t - 2.0) / 2.0;
        _N[ 5] = s * (t + 1.0) * (2.0 * s + t - 2.0) / 2.0;
        _N[ 6] = 2.0 * r * (1.0 - r - s) * (1.0 - t);
        _N[ 7] = 2.0 * r * s * (1.0 - t);
        _N[ 8] = 2.0 * s * (1.0 - r - s) * (1.0 - t);
        _N[ 9] = 2.0 * r * (1.0 - r - s) * (1.0 + t);
        _N[10] = 2.0 * r * s * (1.0 + t);
        _N[11] = 2.0 * s * (1.0 - r - s) * (1.0 + t);
        _N[12] = (1.0 - r - s) * (1.0 - t * t);
        _N[13] = r * (1.0 - t * t);
        _N[14] = s * (1.0 - t * t);

        _dN[0 * nodes +  0] = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
        _dN[0 * nodes +  1] = ((t - 1.0) * (t - 2.0 * r + 2.0)) / 2.0 - r * (t - 1.0);
        _dN[0 * nodes +  2] = 0.0;
        _dN[0 * nodes +  3] = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
        _dN[0 * nodes +  4] = r * (t + 1.0) + ((t + 1.0) * (2.0 * r + t - 2.0)) / 2.0;
        _dN[0 * nodes +  5] = 0.0;
        _dN[0 * nodes +  6] = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * r * (t - 1.0);
        _dN[0 * nodes +  7] = (-2.0) * s * (t - 1.0);
        _dN[0 * nodes +  8] = 2.0 * s * (t - 1.0);
        _dN[0 * nodes +  9] = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * r * (t + 1.0);
        _dN[0 * nodes + 10] =  2.0 * s * (t + 1.0);
        _dN[0 * nodes + 11] =  -2.0 * s * (t + 1.0);
        _dN[0 * nodes + 12] =  t * t - 1.0;
        _dN[0 * nodes + 13] =  1.0 - t * t;
        _dN[0 * nodes + 14] =  0.0;

        _dN[1 * nodes +  0] = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
        _dN[1 * nodes +  1] = 0.0;
        _dN[1 * nodes +  2] = ((t - 1.0) * (t - 2.0 * s + 2.0)) / 2.0 - s * (t - 1.0);
        _dN[1 * nodes +  3] = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
        _dN[1 * nodes +  4] = 0.0;
        _dN[1 * nodes +  5] = s * (t + 1.0) + ((t + 1.0) * (2.0 * s + t - 2.0)) / 2.0;
        _dN[1 * nodes +  6] = 2.0 * r * (t - 1.0);
        _dN[1 * nodes +  7] = (-2.0) * r * (t - 1.0);
        _dN[1 * nodes +  8] = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * s * (t - 1.0);
        _dN[1 * nodes +  9] = (-2.0) * r * (t + 1.0);
        _dN[1 * nodes + 10] =  2.0 * r * (t + 1.0);
        _dN[1 * nodes + 11] =  -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * s * (t + 1.0);
        _dN[1 * nodes + 12] =  t * t - 1.0;
        _dN[1 * nodes + 13] =  0.0;
        _dN[1 * nodes + 14] =  1.0 - t * t;

        _dN[2 * nodes +  0] = -((r + s - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0 - ((t - 1.0) * (r + s - 1.0)) / 2.0;
        _dN[2 * nodes +  1] = (r * (t - 2.0 * r + 2.0)) / 2.0 + (r * (t - 1.0)) / 2.0;
        _dN[2 * nodes +  2] = (s * (t - 2.0 * s + 2.0)) / 2.0 + (s * (t - 1.0)) / 2.0;
        _dN[2 * nodes +  3] = ((r + s - 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 - ((t + 1.0) * (r + s - 1.0)) / 2.0;
        _dN[2 * nodes +  4] = (r * (2.0 * r + t - 2.0)) / 2.0 + (r * (t + 1.0)) / 2.0;
        _dN[2 * nodes +  5] = (s * (2.0 * s + t - 2.0)) / 2.0 + (s * (t + 1.0)) / 2.0;
        _dN[2 * nodes +  6] = 2.0 * r * (r + s - 1.0);
        _dN[2 * nodes +  7] = (-2.0) * r * s;
        _dN[2 * nodes +  8] = 2.0 * s * (r + s - 1.0);
        _dN[2 * nodes +  9] = (-2.0) * r * (r + s - 1.0);
        _dN[2 * nodes + 10] =  2.0 * r * s;
        _dN[2 * nodes + 11] =  (-2.0) * s * (r + s - 1.0);
        _dN[2 * nodes + 12] =  2.0 * t * (r + s - 1.0);
        _dN[2 * nodes + 13] =  (-2.0) * r * t;
        _dN[2 * nodes + 14] =  (-2.0) * s * t;
    }

    static void set()
    {
        double v1 = 1.0 / 6.0;
        double v2 = 4.0 / 6.0;
        double v3 = sqrt(3.0 / 5.0);
        double v4 = 0.0;
        double r[9] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
        double s[9] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
        double t[9] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };

        w[ 0] = 5.0 / 54.0;
        w[ 1] = 5.0 / 54.0;
        w[ 2] = 5.0 / 54.0;
        w[ 3] = 8.0 / 54.0;
        w[ 4] = 8.0 / 54.0;
        w[ 5] = 8.0 / 54.0;
        w[ 6] = 5.0 / 54.0;
        w[ 7] = 5.0 / 54.0;
        w[ 8] = 5.0 / 54.0;
        for (int gp = 0; gp < gps; gp++) {
            set(N + gp * nodes, dN + edim * gp * nodes, r[gp], s[gp], t[gp]);
        }
    }
};

template<>
struct GaussPoints<Element::CODE::HEXA20, 20, 8, 3> {
    constexpr static int nodes = 20, gps = 8, edim = 3;
    static double w[gps], N[gps * nodes], dN[gps * nodes * edim];
    static double cw, cN[nodes], cdN[nodes * edim];

    static void set(double *_N, double *_dN, double r, double s, double t)
    {
        _N[ 0] = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 - t) * (-r - s - t - 2.0));
        _N[ 1] = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 - t) * ( r - s - t - 2.0));
        _N[ 2] = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 - t) * ( r + s - t - 2.0));
        _N[ 3] = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 - t) * (-r + s - t - 2.0));
        _N[ 4] = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 + t) * (-r - s + t - 2.0));
        _N[ 5] = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 + t) * ( r - s + t - 2.0));
        _N[ 6] = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 + t) * ( r + s + t - 2.0));
        _N[ 7] = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 + t) * (-r + s + t - 2.0));
        _N[ 8] = 0.25 * ((1.0 - r * r) * (1.0 - s) * (1.0 - t));
        _N[ 9] = 0.25 * ((1.0 + r) * (1.0 - s * s) * (1.0 - t));
        _N[10] = 0.25 * ((1.0 - r * r) * (1.0 + s) * (1.0 - t));
        _N[11] = 0.25 * ((1.0 - r) * (1.0 - s * s) * (1.0 - t));
        _N[12] = 0.25 * ((1.0 - r * r) * (1.0 - s) * (1.0 + t));
        _N[13] = 0.25 * ((1.0 + r) * (1.0 - s * s) * (1.0 + t));
        _N[14] = 0.25 * ((1.0 - r * r) * (1.0 + s) * (1.0 + t));
        _N[15] = 0.25 * ((1.0 - r) * (1.0 - s * s) * (1.0 + t));
        _N[16] = 0.25 * ((1.0 - r) * (1.0 - s) * (1.0 - t * t));
        _N[17] = 0.25 * ((1.0 + r) * (1.0 - s) * (1.0 - t * t));
        _N[18] = 0.25 * ((1.0 + r) * (1.0 + s) * (1.0 - t * t));
        _N[19] = 0.25 * ((1.0 - r) * (1.0 + s) * (1.0 - t * t));

        _dN[0 * nodes +  0] =  ((s - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
        _dN[0 * nodes +  1] =  ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0 - ((s - 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0;
        _dN[0 * nodes +  2] = -((s + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
        _dN[0 * nodes +  3] = -((s + 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
        _dN[0 * nodes +  4] = -((s - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
        _dN[0 * nodes +  5] = -((s - 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
        _dN[0 * nodes +  6] =  ((s + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
        _dN[0 * nodes +  7] =  ((s + 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 + ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
        _dN[0 * nodes +  8] =  -(r * (s - 1.0) * (t - 1.0)) / 2.0;
        _dN[0 * nodes +  9] =   ((s * s - 1.0) * (t - 1.0)) / 4.0;
        _dN[0 * nodes + 10] =   (r * (s + 1.0) * (t - 1.0)) / 2.0;
        _dN[0 * nodes + 11] =  -((s * s - 1.0) * (t - 1.0)) / 4.0;
        _dN[0 * nodes + 12] =   (r * (s - 1.0) * (t + 1.0)) / 2.0;
        _dN[0 * nodes + 13] =  -((s * s - 1.0) * (t + 1.0)) / 4.0;
        _dN[0 * nodes + 14] =  -(r * (s + 1.0) * (t + 1.0)) / 2.0;
        _dN[0 * nodes + 15] =   ((s * s - 1.0) * (t + 1.0)) / 4.0;
        _dN[0 * nodes + 16] =  -((t * t - 1.0) * (s - 1.0)) / 4.0;
        _dN[0 * nodes + 17] =   ((t * t - 1.0) * (s - 1.0)) / 4.0;
        _dN[0 * nodes + 18] =  -((t * t - 1.0) * (s + 1.0)) / 4.0;
        _dN[0 * nodes + 19] =   ((t * t - 1.0) * (s + 1.0)) / 4.0;

        _dN[1 * nodes +  0] =  ((r - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
        _dN[1 * nodes +  1] = -((r + 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
        _dN[1 * nodes +  2] = -((r + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
        _dN[1 * nodes +  3] =  ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r - 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0;
        _dN[1 * nodes +  4] = -((r - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
        _dN[1 * nodes +  5] =  ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r + 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0;
        _dN[1 * nodes +  6] =  ((r + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
        _dN[1 * nodes +  7] =  ((r - 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
        _dN[1 * nodes +  8] =  -((r * r - 1.0) * (t - 1.0)) / 4.0;
        _dN[1 * nodes +  9] =   (s * (r + 1.0) * (t - 1.0)) / 2.0;
        _dN[1 * nodes + 10] =   ((r * r - 1.0) * (t - 1.0)) / 4.0;
        _dN[1 * nodes + 11] =  -(s * (r - 1.0) * (t - 1.0)) / 2.0;
        _dN[1 * nodes + 12] =   ((r * r - 1.0) * (t + 1.0)) / 4.0;
        _dN[1 * nodes + 13] =  -(s * (r + 1.0) * (t + 1.0)) / 2.0;
        _dN[1 * nodes + 14] =  -((r * r - 1.0) * (t + 1.0)) / 4.0;
        _dN[1 * nodes + 15] =   (s * (r - 1.0) * (t + 1.0)) / 2.0;
        _dN[1 * nodes + 16] =  -((t * t - 1.0) * (r - 1.0)) / 4.0;
        _dN[1 * nodes + 17] =   ((t * t - 1.0) * (r + 1.0)) / 4.0;
        _dN[1 * nodes + 18] =  -((t * t - 1.0) * (r + 1.0)) / 4.0;
        _dN[1 * nodes + 19] =   ((t * t - 1.0) * (r - 1.0)) / 4.0;

        _dN[2 * nodes +  0] =  ((r - 1.0) * (s - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
        _dN[2 * nodes +  1] = -((r + 1.0) * (s - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
        _dN[2 * nodes +  2] =  ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r + 1.0) * (s + 1.0) * (r + s - t - 2.0)) / 8.0;
        _dN[2 * nodes +  3] = -((r - 1.0) * (s + 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
        _dN[2 * nodes +  4] =  ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r - 1.0) * (s - 1.0) * (r + s - t + 2.0)) / 8.0;
        _dN[2 * nodes +  5] = -((r + 1.0) * (s - 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
        _dN[2 * nodes +  6] =  ((r + 1.0) * (s + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
        _dN[2 * nodes +  7] =  ((r - 1.0) * (s + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
        _dN[2 * nodes +  8] =  -((r * r - 1.0) * (s - 1.0)) / 4.0;
        _dN[2 * nodes +  9] =   ((s * s - 1.0) * (r + 1.0)) / 4.0;
        _dN[2 * nodes + 10] =   ((r * r - 1.0) * (s + 1.0)) / 4.0;
        _dN[2 * nodes + 11] =  -((s * s - 1.0) * (r - 1.0)) / 4.0;
        _dN[2 * nodes + 12] =   ((r * r - 1.0) * (s - 1.0)) / 4.0;
        _dN[2 * nodes + 13] =  -((s * s - 1.0) * (r + 1.0)) / 4.0;
        _dN[2 * nodes + 14] =  -((r * r - 1.0) * (s + 1.0)) / 4.0;
        _dN[2 * nodes + 15] =   ((s * s - 1.0) * (r - 1.0)) / 4.0;
        _dN[2 * nodes + 16] =  -(t * (r - 1.0) * (s - 1.0)) / 2.0;
        _dN[2 * nodes + 17] =   (t * (r + 1.0) * (s - 1.0)) / 2.0;
        _dN[2 * nodes + 18] =  -(t * (r + 1.0) * (s + 1.0)) / 2.0;
        _dN[2 * nodes + 19] =   (t * (r - 1.0) * (s + 1.0)) / 2.0;
    }

    static void set()
    {
        double v = 0.577350269189625953;
        double r[8] = {  v,  v,  v,  v, -v, -v, -v, -v };
        double s[8] = { -v, -v,  v,  v, -v, -v,  v,  v };
        double t[8] = { -v,  v, -v,  v, -v,  v, -v,  v };

        for (int gp = 0; gp < gps; gp++) {
            w[gp] = 1;
            set(N + gp * nodes, dN + edim * gp * nodes, r[gp], s[gp], t[gp]);
        }
    }
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_BASIS_H_ */
