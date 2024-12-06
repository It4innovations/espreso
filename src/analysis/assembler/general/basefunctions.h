
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_BASEFUNCTIONS_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_BASEFUNCTIONS_H_

#include "element.h"
#include "mesh/element.h"

#include <array>

namespace espreso {

template <Element::CODE code, size_t gps> struct BaseFunctions;

template<>
struct BaseFunctions<Element::CODE::TRIANGLE3, 6> {
    constexpr static size_t nodes = 3, gps = 6, edim = 2;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s)
    {
         N[0]    = 1 - r - s;
         N[1]    = r;
         N[2]    = s;

        dN[0][0] = -1;
        dN[1][0] =  1;
        dN[2][0] =  0;

        dN[0][1] = -1;
        dN[1][1] =  0;
        dN[2][1] =  1;
    }

    template <typename Element>
    static void simd(Element &element)
    {
        double w[gps] = { 0.111690794839005, 0.111690794839005, 0.111690794839005, 0.054975871827661, 0.054975871827661, 0.054975871827661 };
        double r[gps] = { 0.445948490915965, 0.445948490915965, 0.108103018168070, 0.091576213509771, 0.091576213509771, 0.816847572980459 };
        double s[gps] = { 0.445948490915965, 0.108103018168070, 0.445948490915965, 0.091576213509771, 0.816847572980459, 0.091576213509771 };

        for (size_t gp = 0; gp < gps; gp++) {
            element.w[gp] = w[gp];
            set(element.N[gp], element.dN[gp], r[gp], s[gp]);
        }

        double nn[2] = { 0, 1 };
        set(element.NN[0], element.dNN[0], nn[0], nn[0]);
        set(element.NN[1], element.dNN[1], nn[0], nn[1]);
        set(element.NN[2], element.dNN[2], nn[1], nn[0]);
    }
};

template<>
struct BaseFunctions<Element::CODE::TRIANGLE6, 6> {
    constexpr static size_t nodes = 6, gps = 6, edim = 2;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s)
    {
         N[0]    = (1.0 - r - s) * (1.0 - 2.0 * (r + s));
         N[1]    = -(r) * (1.0 - 2.0 * r);
         N[2]    = -(s) * (1.0 - 2.0 * s);
         N[3]    = 4.0 * (r) * (1.0 - r - s);
         N[4]    = 4.0 * (r) * (s);
         N[5]    = 4.0 * (s) * (1.0 - r - s);

        dN[0][0] = -3.0 + 4.0 * r + 4.0 * s;
        dN[1][0] = -1.0 + 4.0 * r;
        dN[2][0] = 0.0;
        dN[3][0] = 4.0 - 8.0 * r - 4.0 * s;
        dN[4][0] = 4.0 * s;
        dN[5][0] = -4.0 * s;

        dN[0][1] = -3.0 + 4.0 * r + 4.0 * s;
        dN[1][1] = 0.0;
        dN[2][1] = -1.0 + 4.0 * s;
        dN[3][1] = -4.0 * r;
        dN[4][1] = 4.0 * r;
        dN[5][1] = 4.0 - 4.0 * r - 8.0 * s;
    }

    template <typename Element>
    static void simd(Element &element)
    {
        double w[gps] = { 0.109951743655322, 0.109951743655322, 0.109951743655322, 0.223381589678011, 0.223381589678011, 0.223381589678011 };
        double r[gps] = { 0.091576213509771, 0.816847572980459, 0.091576213509771, 0.445948490915965, 0.108103018168070, 0.445948490915965 };
        double s[gps] = { 0.091576213509771, 0.091576213509771, 0.816847572980459, 0.445948490915965, 0.445948490915965, 0.108103018168070 };

        for (size_t gp = 0; gp < gps; gp++) {
            element.w[gp] = .5 * w[gp];
            set(element.N[gp], element.dN[gp], r[gp], s[gp]);
        }

        double nn[3] = { 0, 1, .5 };
        set(element.NN[0], element.dNN[0], nn[0], nn[0]);
        set(element.NN[1], element.dNN[1], nn[1], nn[0]);
        set(element.NN[2], element.dNN[2], nn[0], nn[1]);
        set(element.NN[3], element.dNN[3], nn[2], nn[0]);
        set(element.NN[4], element.dNN[4], nn[2], nn[2]);
        set(element.NN[5], element.dNN[5], nn[0], nn[2]);
    }
};

template<>
struct BaseFunctions<Element::CODE::TETRA4, 4> {
    constexpr static size_t nodes = 4, gps = 4, edim = 3;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s, double t)
    {
         N[0]    = r;
         N[1]    = s;
         N[2]    = t;
         N[3]    = 1.0 - r - s - t;

        dN[0][0] =  1.0;
        dN[1][0] =  0.0;
        dN[2][0] =  0.0;
        dN[3][0] = -1.0;

        dN[0][1] =  0.0;
        dN[1][1] =  0.0;
        dN[2][1] =  1.0;
        dN[3][1] = -1.0;

        dN[0][2] =  0.0;
        dN[1][2] =  1.0;
        dN[2][2] =  0.0;
        dN[3][2] = -1.0;
    }

    template <typename Element>
    static void simd(Element &element)
    {
        double r[gps] = { 0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105 };
        double s[gps] = { 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685 };
        double t[gps] = { 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105 };

        for (size_t gp = 0; gp < gps; ++gp) {
            element.w[gp] = 1.0 / 24.0;
            set(element.N[gp], element.dN[gp], r[gp], s[gp], t[gp]);
        }

        double nn[2] = { 0, 1 };
        set(element.NN[0], element.dNN[0], nn[1], nn[0], nn[0]);
        set(element.NN[1], element.dNN[1], nn[0], nn[1], nn[0]);
        set(element.NN[2], element.dNN[2], nn[0], nn[0], nn[1]);
        set(element.NN[3], element.dNN[3], nn[0], nn[0], nn[0]);
    }
};

template<>
struct BaseFunctions<Element::CODE::TETRA10, 15> {
    constexpr static size_t nodes = 10, gps = 15, edim = 3;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s, double t)
    {
         N[0]    = r * (2.0 * r - 1.0);
         N[1]    = s * (2.0 * s - 1.0);
         N[2]    = t * (2.0 * t - 1.0);
         N[3]    = 2.0 * r * r + 4.0 * r * s + 4.0 * r * t - 3.0 * r + 2.0* s * s + 4.0 * s * t - 3.0 * s + 2.0 * t * t - 3.0 * t + 1.0;
         N[4]    = 4.0 * r * s;
         N[5]    = 4.0 * s * t;
         N[6]    = 4.0 * r * t;
         N[7]    = r * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
         N[8]    = s * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);
         N[9]    = t * (-4.0 * r - 4.0 * s - 4.0 * t + 4.0);

        dN[0][0] = 4.0 * r - 1.0;
        dN[1][0] = 0;
        dN[2][0] = 0;
        dN[3][0] = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
        dN[4][0] = 4.0 * s;
        dN[5][0] = 0;
        dN[6][0] = 4.0 * t;
        dN[7][0] = -8.0 * r - 4.0 * s - 4.0 * t + 4.0;
        dN[8][0] = -4.0 * s;
        dN[9][0] = -4.0 * t;

        dN[0][1] = 0;
        dN[1][1] = 0;
        dN[2][1] = 4.0 * t - 1.0;
        dN[3][1] = 4.0 * r + 4.0 * s + 4.0* t  - 3.0;
        dN[4][1] = 0;
        dN[5][1] = 4.0 * s;
        dN[6][1] = 4.0 * r;
        dN[7][1] = -4.0 * r;
        dN[8][1] = -4.0 * s;
        dN[9][1] = -4.0 * r - 4.0 * s - 8.0 * t + 4.0;

        dN[0][2] = 0;
        dN[1][2] = 4.0 * s - 1.0;
        dN[2][2] = 0 ;
        dN[3][2] = 4.0 * r + 4.0 * s + 4.0 * t - 3.0;
        dN[4][2] = 4.0 * r;
        dN[5][2] = 4.0 * t;
        dN[6][2] = 0;
        dN[7][2] = -4.0 * r;
        dN[8][2] = -4.0 * r - 8.0 * s - 4.0 * t + 4.0;
        dN[9][2] = -4.0 * t;
    }

    template <typename Element>
    static void simd(Element &element)
    {
        double w[gps] = {
                0.030283678097089 , 0.006026785714286 , 0.006026785714286 , 0.006026785714286 ,
                0.006026785714286 , 0.011645249086029 , 0.011645249086029 , 0.011645249086029 ,
                0.011645249086029 , 0.010949141561386 , 0.010949141561386 , 0.010949141561386 ,
                0.010949141561386 , 0.010949141561386 , 0.010949141561386  };
        double r[gps] = {
                0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333,
                0.3333333333333333, 0.7272727272727273, 0.0909090909090909, 0.0909090909090909,
                0.0909090909090909, 0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
                0.0665501535736643, 0.4334498464263357, 0.4334498464263357 };
        double s[gps] = {
                0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333,
                0.0000000000000000, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909,
                0.7272727272727273, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
                0.4334498464263357, 0.0665501535736643, 0.4334498464263357 };
        double t[gps] = {
                0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000,
                0.3333333333333333, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273,
                0.0909090909090909, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
                0.4334498464263357, 0.4334498464263357, 0.0665501535736643 };

        for (size_t gp = 0; gp < gps; ++gp) {
            element.w[gp] = w[gp];
            set(element.N[gp], element.dN[gp], r[gp], s[gp], t[gp]);
        }

        double nn[3] = { 0, 1, .5 };
        set(element.NN[0], element.dNN[0], nn[1], nn[0], nn[0]);
        set(element.NN[1], element.dNN[1], nn[0], nn[1], nn[0]);
        set(element.NN[2], element.dNN[2], nn[0], nn[0], nn[1]);
        set(element.NN[3], element.dNN[3], nn[0], nn[0], nn[0]);

        set(element.NN[4], element.dNN[4], nn[2], nn[2], nn[0]);
        set(element.NN[5], element.dNN[5], nn[0], nn[2], nn[2]);
        set(element.NN[6], element.dNN[6], nn[2], nn[0], nn[2]);
        set(element.NN[7], element.dNN[7], nn[2], nn[0], nn[0]);
        set(element.NN[8], element.dNN[8], nn[0], nn[2], nn[0]);
        set(element.NN[9], element.dNN[9], nn[0], nn[0], nn[2]);
    }
};

template<>
struct BaseFunctions<Element::CODE::PYRAMID13, 14> {
    constexpr static size_t nodes = 13, gps = 14, edim = 3;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s, double t)
    {
         N[ 0]    = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 - s) * (-1.0 - (0.5 * (1.0 - t)) * r - (0.5 * (1.0 - t)) * s));
         N[ 1]    = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 - s) * (-1.0 + (0.5 * (1 - t))   * r - (0.5 * (1.0 - t)) * s));
         N[ 2]    = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 + s) * (-1.0 + (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
         N[ 3]    = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 + s) * (-1.0 - (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
         N[ 4]    = (1.0 - (0.5 * (1.0 - t))) * (1.0 - 2.0 * (0.5 * (1.0 - t)));
         N[ 5]    = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - s) * (1.0 - r * r);
         N[ 6]    = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + r) * (1.0 - s * s);
         N[ 7]    = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + s) * (1.0 - r * r);
         N[ 8]    = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - r) * (1.0 - s * s);
         N[ 9]    = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r - s + r * s);
         N[10]    = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r - s - r * s);
         N[11]    = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r + s + r * s);
         N[12]    = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r + s - r * s);

        dN[ 0][0] = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
        dN[ 1][0] = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
        dN[ 2][0] =  (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0) + (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
        dN[ 3][0] =  (t / 2.0 - 0.5) * (s + 1.0) * (r - 1.0) * (t / 8.0 - 1.0 / 8.0) - (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0);
        dN[ 4][0] =  0.0;
        dN[ 5][0] =  r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s - 1.0);
        dN[ 6][0] = -((s * s - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
        dN[ 7][0] = -r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s + 1.0);
        dN[ 8][0] =  ((s * s - 1) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
        dN[ 9][0] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
        dN[10][0] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
        dN[11][0] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);
        dN[12][0] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);

        dN[ 0][1] = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) - 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
        dN[ 1][1] =  (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0) - (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) - s * (t / 2.0 - 1.0 / 2.0) + 1.0);
        dN[ 2][1] =  (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) + 1.0) + (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
        dN[ 3][1] = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s * (t / 2.0 - 1.0 / 2.0) - r * (t / 2.0 - 1.0 / 2.0) + 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0);
        dN[ 4][1] =  0.0;
        dN[ 5][1] =  ((r * r - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
        dN[ 6][1] = -s * (t / 2.0 - 1.0 / 2.0) * (t / 2.0 - 0.5) * (r + 1.0);
        dN[ 7][1] = -((r * r - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
        dN[ 8][1] =  s * (t / 2.0 - 0.5) * (t / 2.0 - 0.5) * (r - 1.0);
        dN[ 9][1] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);
        dN[10][1] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
        dN[11][1] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
        dN[12][1] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);

        dN[ 0][2] = -((r - 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0)) / 8.0 - (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
        dN[ 1][2] = -((r + 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0)) / 8.0 - (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
        dN[ 2][2] =  ((r + 1.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0)) / 8.0 + (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
        dN[ 3][2] =  (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0) - ((r - 1.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0)) / 8.0;
        dN[ 4][2] =  t + 0.5;
        dN[ 5][2] =  ((r * r - 1.0) * (t / 2.0 - 0.5) * (s - 1.0)) / 2.0;
        dN[ 6][2] = -((s * s - 1.0) * (t / 2.0 - 0.5) * (r + 1.0)) / 2.0;
        dN[ 7][2] = -((r * r - 1.0) * (t / 2.0 - 0.5) * (s + 1.0)) / 2.0;
        dN[ 8][2] =  ((s * s - 1.0) * (t / 2.0 - 0.5) * (r - 1.0)) / 2.0;
        dN[ 9][2] =  ((t / 2.0 - 0.5) * (r + s - r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r + s - r * s - 1.0)) / 2.0;
        dN[10][2] = -((t / 2.0 - 0.5) * (r - s - r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r - s - r * s + 1.0)) / 2.0;
        dN[11][2] = -((t / 2.0 - 0.5) * (r + s + r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r + s + r * s + 1.0)) / 2.0;
        dN[12][2] =  ((t / 2.0 - 0.5) * (r - s + r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r - s + r * s - 1.0)) / 2.0;
    }

    template <typename Element>
    static void simd(Element &element)
    {
        double v1 = 0.758786910639329015;
        double v2 = 0.795822425754222018;
        double v3 = 0;
        double w1 = 0.335180055401662;
        double w2 = 0.886426592797784;
        double w[gps] = {  w1,  w1,  w1,  w1,  w1,  w1,  w1,  w1,  w2,  w2,  w2, w2,  w2, w2 };
        double r[gps] = { -v1,  v1,  v1, -v1, -v1,  v1,  v1, -v1,  v3,  v3,  v2, v3, -v2, v3 };
        double s[gps] = { -v1, -v1,  v1,  v1, -v1, -v1,  v1,  v1,  v3, -v2,  v3, v2,  v3, v3 };
        double t[gps] = { -v1, -v1, -v1, -v1,  v1,  v1,  v1,  v1, -v2,  v3,  v3, v3,  v3, v2 };

        for (size_t gp = 0; gp < gps; gp++) {
            element.w[gp] = w[gp];
            set(element.N[gp], element.dN[gp], r[gp], s[gp], t[gp]);
        }

        double nn[3] = { -1, 1, 0 };
        set(element.NN[ 0], element.dNN[ 0], nn[0], nn[0], nn[0]);
        set(element.NN[ 1], element.dNN[ 1], nn[1], nn[0], nn[0]);
        set(element.NN[ 2], element.dNN[ 2], nn[1], nn[1], nn[0]);
        set(element.NN[ 3], element.dNN[ 3], nn[0], nn[1], nn[0]);
        set(element.NN[ 4], element.dNN[ 4], nn[2], nn[2], nn[1]);

        set(element.NN[ 5], element.dNN[ 5], nn[2], nn[0], nn[0]);
        set(element.NN[ 6], element.dNN[ 6], nn[1], nn[2], nn[0]);
        set(element.NN[ 7], element.dNN[ 7], nn[2], nn[1], nn[0]);
        set(element.NN[ 8], element.dNN[ 8], nn[0], nn[2], nn[0]);
        set(element.NN[ 9], element.dNN[ 9], nn[0], nn[0], nn[2]);
        set(element.NN[10], element.dNN[10], nn[1], nn[0], nn[2]);
        set(element.NN[11], element.dNN[11], nn[1], nn[1], nn[2]);
        set(element.NN[12], element.dNN[12], nn[0], nn[1], nn[2]);
    }
};

template<>
struct BaseFunctions<Element::CODE::PRISMA6, 9> {
    constexpr static size_t nodes = 6, gps = 9, edim = 3;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s, double t)
    {
         N[0]    = 0.5 * ((1.0 - t) * (1.0 - r - s));
         N[1]    = 0.5 * ((1.0 - t) * r);
         N[2]    = 0.5 * ((1.0 - t) * s);
         N[3]    = 0.5 * ((1.0 + t) * (1.0 - r - s));
         N[4]    = 0.5 * ((1.0 + t) * r);
         N[5]    = 0.5 * ((1.0 + t) * s);

        dN[0][0] =  t / 2.0 - 1.0 / 2.0;
        dN[1][0] = -t / 2.0 + 1.0 / 2.0;
        dN[2][0] =  0.0;
        dN[3][0] = -t / 2.0 - 1.0 / 2.0;
        dN[4][0] =  t / 2.0 + 1.0 / 2.0;
        dN[5][0] =  0;

        dN[0][1] =  t / 2.0 - 1.0 / 2.0;
        dN[1][1] =  0.0;
        dN[2][1] = -t / 2.0 + 1.0 / 2.0;
        dN[3][1] = -t / 2.0 - 1.0 / 2.0;
        dN[4][1] =  0.0;
        dN[5][1] =  t / 2.0 + 1.0 / 2.0;

        dN[0][2] =  r / 2.0 + s / 2.0 - 1.0 / 2.0;
        dN[1][2] = -r / 2.0;
        dN[2][2] =          - s / 2.0;
        dN[3][2] = -r / 2.0 - s / 2.0 + 1.0 / 2.0;
        dN[4][2] =  r / 2.0;
        dN[5][2] =            s / 2.0;
    }

    template <typename Element>
    static void simd(Element &element)
    {
        double v1 = 1.0 / 6.0;
        double v2 = 4.0 / 6.0;
        double v3 = sqrt(3.0 / 5.0);
        double v4 = 0.0;
        double w1 = 5.0 / 54.;
        double w2 = 8.0 / 54.;
        double w[gps] = {  w1,  w1,  w1,  w2,  w2,  w2,  w1,  w1,  w1 };
        double r[gps] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
        double s[gps] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
        double t[gps] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };

        for (size_t gp = 0; gp < gps; gp++) {
            element.w[gp] = w[gp];
            set(element.N[gp], element.dN[gp], r[gp], s[gp], t[gp]);
        }

        double nn[3] = { -1, 0, 1 };
        set(element.NN[0], element.dNN[0], nn[1], nn[1], nn[0]);
        set(element.NN[1], element.dNN[1], nn[2], nn[1], nn[0]);
        set(element.NN[2], element.dNN[2], nn[1], nn[2], nn[0]);
        set(element.NN[3], element.dNN[3], nn[1], nn[1], nn[2]);
        set(element.NN[4], element.dNN[4], nn[2], nn[1], nn[2]);
        set(element.NN[5], element.dNN[5], nn[1], nn[2], nn[2]);
    }
};

template<>
struct BaseFunctions<Element::CODE::PRISMA15, 9> {
    constexpr static size_t nodes = 15, gps = 9, edim = 3;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s, double t)
    {
         N[ 0]    = -(1.0 - r - s) * (1.0 - t) * (2.0 * r + 2.0 * s + t) / 2.0;
         N[ 1]    = r * (1.0 - t) * (2.0 * r - t - 2.0) / 2.0;
         N[ 2]    = s * (1.0 - t) * (2.0 * s - t - 2.0) / 2.0;
         N[ 3]    = -(1.0 - r - s) * (1.0 + t) * (2.0 * r + 2.0 * s - t) / 2.0;
         N[ 4]    = r * (t + 1.0) * (2.0 * r + t - 2.0) / 2.0;
         N[ 5]    = s * (t + 1.0) * (2.0 * s + t - 2.0) / 2.0;
         N[ 6]    = 2.0 * r * (1.0 - r - s) * (1.0 - t);
         N[ 7]    = 2.0 * r * s * (1.0 - t);
         N[ 8]    = 2.0 * s * (1.0 - r - s) * (1.0 - t);
         N[ 9]    = 2.0 * r * (1.0 - r - s) * (1.0 + t);
         N[10]    = 2.0 * r * s * (1.0 + t);
         N[11]    = 2.0 * s * (1.0 - r - s) * (1.0 + t);
         N[12]    = (1.0 - r - s) * (1.0 - t * t);
         N[13]    = r * (1.0 - t * t);
         N[14]    = s * (1.0 - t * t);

        dN[ 0][0] = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
        dN[ 1][0] = ((t - 1.0) * (t - 2.0 * r + 2.0)) / 2.0 - r * (t - 1.0);
        dN[ 2][0] = 0.0;
        dN[ 3][0] = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
        dN[ 4][0] = r * (t + 1.0) + ((t + 1.0) * (2.0 * r + t - 2.0)) / 2.0;
        dN[ 5][0] = 0.0;
        dN[ 6][0] = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * r * (t - 1.0);
        dN[ 7][0] = (-2.0) * s * (t - 1.0);
        dN[ 8][0] = 2.0 * s * (t - 1.0);
        dN[ 9][0] = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * r * (t + 1.0);
        dN[10][0] =  2.0 * s * (t + 1.0);
        dN[11][0] =  -2.0 * s * (t + 1.0);
        dN[12][0] =  t * t - 1.0;
        dN[13][0] =  1.0 - t * t;
        dN[14][0] =  0.0;

        dN[ 0][1] = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
        dN[ 1][1] = 0.0;
        dN[ 2][1] = ((t - 1.0) * (t - 2.0 * s + 2.0)) / 2.0 - s * (t - 1.0);
        dN[ 3][1] = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
        dN[ 4][1] = 0.0;
        dN[ 5][1] = s * (t + 1.0) + ((t + 1.0) * (2.0 * s + t - 2.0)) / 2.0;
        dN[ 6][1] = 2.0 * r * (t - 1.0);
        dN[ 7][1] = (-2.0) * r * (t - 1.0);
        dN[ 8][1] = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * s * (t - 1.0);
        dN[ 9][1] = (-2.0) * r * (t + 1.0);
        dN[10][1] =  2.0 * r * (t + 1.0);
        dN[11][1] =  -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * s * (t + 1.0);
        dN[12][1] =  t * t - 1.0;
        dN[13][1] =  0.0;
        dN[14][1] =  1.0 - t * t;

        dN[ 0][2] = -((r + s - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0 - ((t - 1.0) * (r + s - 1.0)) / 2.0;
        dN[ 1][2] = (r * (t - 2.0 * r + 2.0)) / 2.0 + (r * (t - 1.0)) / 2.0;
        dN[ 2][2] = (s * (t - 2.0 * s + 2.0)) / 2.0 + (s * (t - 1.0)) / 2.0;
        dN[ 3][2] = ((r + s - 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 - ((t + 1.0) * (r + s - 1.0)) / 2.0;
        dN[ 4][2] = (r * (2.0 * r + t - 2.0)) / 2.0 + (r * (t + 1.0)) / 2.0;
        dN[ 5][2] = (s * (2.0 * s + t - 2.0)) / 2.0 + (s * (t + 1.0)) / 2.0;
        dN[ 6][2] = 2.0 * r * (r + s - 1.0);
        dN[ 7][2] = (-2.0) * r * s;
        dN[ 8][2] = 2.0 * s * (r + s - 1.0);
        dN[ 9][2] = (-2.0) * r * (r + s - 1.0);
        dN[10][2] =  2.0 * r * s;
        dN[11][2] =  (-2.0) * s * (r + s - 1.0);
        dN[12][2] =  2.0 * t * (r + s - 1.0);
        dN[13][2] =  (-2.0) * r * t;
        dN[14][2] =  (-2.0) * s * t;
    }

    template <typename Element>
    static void simd(Element &element)
    {
        double v1 = 1.0 / 6.0;
        double v2 = 4.0 / 6.0;
        double v3 = sqrt(3.0 / 5.0);
        double v4 = 0.0;
        double w1 = 5.0 / 54.0;
        double w2 = 8.0 / 54.0;
        double w[gps] = {  w1,  w1,  w1,  w2,  w2,  w2,  w1,  w1,  w1 };
        double r[gps] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
        double s[gps] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
        double t[gps] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };

        for (size_t gp = 0; gp < gps; gp++) {
            element.w[gp] = w[gp];
            set(element.N[gp], element.dN[gp], r[gp], s[gp], t[gp]);
        }

        double nn[4] = { -1, 0, 1, .5 };
        set(element.NN[ 0], element.dNN[ 0], nn[1], nn[1], nn[0]);
        set(element.NN[ 1], element.dNN[ 1], nn[2], nn[1], nn[0]);
        set(element.NN[ 2], element.dNN[ 2], nn[1], nn[2], nn[0]);
        set(element.NN[ 3], element.dNN[ 3], nn[1], nn[1], nn[2]);
        set(element.NN[ 4], element.dNN[ 4], nn[2], nn[1], nn[2]);
        set(element.NN[ 5], element.dNN[ 5], nn[1], nn[2], nn[2]);

        set(element.NN[ 6], element.dNN[ 6], nn[3], nn[1], nn[0]);
        set(element.NN[ 7], element.dNN[ 7], nn[3], nn[3], nn[0]);
        set(element.NN[ 8], element.dNN[ 8], nn[1], nn[3], nn[0]);
        set(element.NN[ 9], element.dNN[ 9], nn[3], nn[1], nn[2]);
        set(element.NN[10], element.dNN[10], nn[3], nn[3], nn[2]);
        set(element.NN[11], element.dNN[11], nn[1], nn[3], nn[2]);
        set(element.NN[12], element.dNN[12], nn[1], nn[1], nn[1]);
        set(element.NN[13], element.dNN[13], nn[2], nn[1], nn[1]);
        set(element.NN[14], element.dNN[14], nn[1], nn[2], nn[1]);
    }
};

double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
{
    return curr == prev ? curr : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}

double constexpr sqrt(double x)
{
    return sqrtNewtonRaphson(x, x, 0);
}

template <size_t gps> struct GGaussePoints;

template <> struct GGaussePoints<1> {

    constexpr static std::array<double, 1> gps{ 0 };
    constexpr static std::array<double, 1> w  { 2 };
};

template <> struct GGaussePoints<2> {
    constexpr static double qa = 1 / sqrt(3);

    constexpr static std::array<double, 2> gps{ -qa, qa };
    constexpr static std::array<double, 2> w  {   1,  1 };
};

template <> struct GGaussePoints<3> {
    constexpr static double wa = 5 / 9.;
    constexpr static double wb = 8 / 9.;
    constexpr static double qa = sqrt(3 / 5.);
    constexpr static double qb = 0;

    constexpr static std::array<double, 3> gps{ -qa, qb, qa };
    constexpr static std::array<double, 3> w  {  wa, wb, wa };
};

template <> struct GGaussePoints<4> {
    constexpr static double wa = (18. - sqrt(30)) / 36.;
    constexpr static double wb = (18. + sqrt(30)) / 36.;
    constexpr static double qa = sqrt(3 / 7. + sqrt(6 / 5.) * 2 / 7.);
    constexpr static double qb = sqrt(3 / 7. - sqrt(6 / 5.) * 2 / 7.);

    constexpr static std::array<double, 4> gps{ -qa, -qb, qb, qa };
    constexpr static std::array<double, 4> w  {  wa,  wb, wb, wa };
};

template <> struct GGaussePoints<5> {
    constexpr static double wa = (322. - 13. * sqrt(70.)) / 900.;
    constexpr static double wb = (322. + 13. * sqrt(70.)) / 900.;
    constexpr static double wc = 128 / 225.;
    constexpr static double qa = sqrt(5. + 2. * sqrt(10 / 7.)) / 3.;
    constexpr static double qb = sqrt(5. - 2. * sqrt(10 / 7.)) / 3.;
    constexpr static double qc = 0;

    constexpr static std::array<double, 5> gps{ -qa, -qb, qc, qb, qa };
    constexpr static std::array<double, 5> w  {  wa,  wb, wc, wb, wa };
};

template <> struct GGaussePoints<6> {
    constexpr static double wa = 0.171324492379170;
    constexpr static double wb = 0.360761573048139;
    constexpr static double wc = 0.467913934572691;
    constexpr static double qa = 0.932469514203152;
    constexpr static double qb = 0.661209386466265;
    constexpr static double qc = 0.238619186083197;

    constexpr static std::array<double, 6> gps{ -qa, -qb, -qc, qc, qb, qa };
    constexpr static std::array<double, 6> w  {  wa,  wb,  wc, wc, wb, wa };
};

template <> struct GGaussePoints<7> {
    constexpr static double wa = 0.129484966168870;
    constexpr static double wb = 0.279705391489277;
    constexpr static double wc = 0.381830050505119;
    constexpr static double wd = 0.417959183673469;
    constexpr static double qa = 0.949107912342759;
    constexpr static double qb = 0.741531185599394;
    constexpr static double qc = 0.405845151377397;
    constexpr static double qd = 0;

    constexpr static std::array<double, 7> gps{ -qa, -qb, -qc, qd, qc, qb, qa };
    constexpr static std::array<double, 7> w  {  wa,  wb,  wc, wd, wc, wb, wa };
};

template <> struct GGaussePoints<8> {
    constexpr static double wa = 0.101228536290376;
    constexpr static double wb = 0.222381034453374;
    constexpr static double wc = 0.313706645877887;
    constexpr static double wd = 0.362683783378362;
    constexpr static double qa = 0.960289856497536;
    constexpr static double qb = 0.796666477413627;
    constexpr static double qc = 0.525532409916329;
    constexpr static double qd = 0.183434642495650;

    constexpr static std::array<double, 8> gps{ -qa, -qb, -qc, -qd, qd, qc, qb, qa };
    constexpr static std::array<double, 8> w  {  wa,  wb,  wc,  wd, wd, wc, wb, wa };
};

template <size_t gps, size_t edim>
constexpr size_t getorder()
{
    switch (edim) {
    case 3:
        switch (gps) {
        case   1: return 1;
        case   8: return 2;
        case  27: return 3;
        case  64: return 4;
        case 125: return 5;
        case 216: return 6;
        case 343: return 7;
        case 512: return 8;
        default: return -1;
        }; break;
    case 2:
        switch (gps) {
        case   1: return 1;
        case   4: return 2;
        case   9: return 3;
        case  16: return 4;
        case  25: return 5;
        case  36: return 6;
        case  49: return 7;
        case  64: return 8;
        default: return -1;
        }; break;
    case 1:
        return gps;
    }
    return -1; // error
}

template<typename E, size_t nodes, size_t edim> struct GaussPointsDegradedSetter;
template<typename E, size_t nodes, size_t edim> struct GaussPointsDegraded;

template <size_t ngps>
struct BaseFunctions<Element::CODE::LINE2 , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::LINE2, ngps>, 2, 1> {
    constexpr static size_t nodes = 2, gps = ngps, edim = 1, order = getorder<gps, edim>();
    constexpr static std::array<int, 2> n_order = { 0, 1 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::LINE3, ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::LINE3, ngps>, 3, 1> {
    constexpr static size_t nodes = 3, gps = ngps, edim = 1, order = getorder<gps, edim>();
    constexpr static std::array<int, 3> n_order = { 0, 1, 2 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::TRIANGLE3, ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::TRIANGLE3, ngps>, 4, 2> {
    constexpr static size_t nodes = 3, gps = ngps, edim = 2, order = getorder<gps, edim>();
    constexpr static std::array<int, 4> n_order = { 0, 1, 2, 2 };
};


template <size_t ngps>
struct BaseFunctions<Element::CODE::SQUARE4  , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::SQUARE4  , ngps>, 4, 2> {
    constexpr static size_t nodes = 4, gps = ngps, edim = 2, order = getorder<gps, edim>();
    constexpr static std::array<int, 4> n_order = { 0, 1, 2, 3 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::TRIANGLE6, ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::TRIANGLE6, ngps>, 8, 2> {
    constexpr static size_t nodes = 6, gps = ngps, edim = 2, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 2, 3, 4, 2, 5 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::SQUARE8  , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::SQUARE8  , ngps>, 8, 2> {
    constexpr static size_t nodes = 8, gps = ngps, edim = 2, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 3, 4, 5, 6, 7 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::TETRA4  , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::TETRA4   , ngps>, 8, 3> {
    constexpr static size_t nodes = 4, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 2, 3, 3, 3, 3 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::PYRAMID5, ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::PYRAMID5 , ngps>, 8, 3> {
    constexpr static size_t nodes = 5, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 3, 4, 4, 4, 4 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::PRISMA6 , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::PRISMA6  , ngps>, 8, 3> {
    constexpr static size_t nodes = 6, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 2, 3, 4, 5, 5 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::HEXA8   , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::HEXA8    , ngps>, 8, 3> {
    constexpr static size_t nodes = 8, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 3, 4, 5, 6, 7 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::TETRA10  , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::TETRA10  , ngps>, 20, 3> {
    constexpr static size_t nodes = 10, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 2, 3, 3, 3, 3, 4, 5, 2, 6, 3, 3, 3, 3, 7, 8, 9, 9 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::PYRAMID13, ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::PYRAMID13, ngps>, 20, 3> {
    constexpr static size_t nodes = 13, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, 4, 4, 4, 4, 9, 10, 11, 12 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::PRISMA15 , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::PRISMA15 , ngps>, 20, 3> {
    constexpr static size_t nodes = 15, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 2, 3, 4, 5, 5, 6, 7, 2, 8, 9, 10, 5, 11, 12, 13, 14, 14 };
};

template <size_t ngps>
struct BaseFunctions<Element::CODE::HEXA20   , ngps>: GaussPointsDegraded<BaseFunctions<Element::CODE::HEXA20   , ngps>, 20, 3> {
    constexpr static size_t nodes = 20, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
};

template<typename E> struct GaussPointsDegradedSetter<E, 2, 1>
{
    static void set(double N[], double dN[][1], double r)
    {
        for (size_t n = 0; n < E::nodes; ++n) {
            N[n] = dN[n][0] = 0;
        }

         N[E::n_order[0]]    += 0.5 * (1 - r);
         N[E::n_order[1]]    += 0.5 * (1 + r);

        dN[E::n_order[0]][0] += 0.5 * (  - 1);
        dN[E::n_order[1]][0] += 0.5 * (  + 1);
    }
};

template<typename E> struct GaussPointsDegradedSetter<E, 3, 1>
{
    static void set(double N[], double dN[][1], double r)
    {
        for (size_t n = 0; n < E::nodes; ++n) {
            N[n] = dN[n][0] = 0;
        }

         N[E::n_order[0]]    += 0.5 * (1 - r) * (-r);
         N[E::n_order[1]]    += 0.5 * (1 + r) * ( r);
         N[E::n_order[2]]    += 1 - r * r;

        dN[E::n_order[0]][0] += 0.5 * (  - 1) * (-r) - 0.5 * (1 - r);
        dN[E::n_order[1]][0] += 0.5 * (  + 1) * ( r) + 0.5 * (1 + r);
        dN[E::n_order[2]][0] +=   - 2 * r;
    }
};

template<typename E> struct GaussPointsDegradedSetter<E, 4, 2>
{
    static void set(double N[], double dN[][2], double r, double s)
    {
        for (size_t n = 0; n < E::nodes; ++n) {
            N[n] = dN[n][0] = dN[n][1] = 0;
        }

         N[E::n_order[0]]    += 0.25 * (1 - r) * (1 - s);
         N[E::n_order[1]]    += 0.25 * (1 + r) * (1 - s);
         N[E::n_order[2]]    += 0.25 * (1 + r) * (1 + s);
         N[E::n_order[3]]    += 0.25 * (1 - r) * (1 + s);

        dN[E::n_order[0]][0] += 0.25 * (  - 1) * (1 - s);
        dN[E::n_order[1]][0] += 0.25 * (  + 1) * (1 - s);
        dN[E::n_order[2]][0] += 0.25 * (  + 1) * (1 + s);
        dN[E::n_order[3]][0] += 0.25 * (  - 1) * (1 + s);

        dN[E::n_order[0]][1] += 0.25 * (1 - r) * (  - 1);
        dN[E::n_order[1]][1] += 0.25 * (1 + r) * (  - 1);
        dN[E::n_order[2]][1] += 0.25 * (1 + r) * (  + 1);
        dN[E::n_order[3]][1] += 0.25 * (1 - r) * (  + 1);
    }
};

template<typename E> struct GaussPointsDegradedSetter<E, 8, 2>
{
    static void set(double N[], double dN[][2], double r, double s)
    {
        for (size_t n = 0; n < E::nodes; ++n) {
            N[n] = dN[n][0] = dN[n][1] = 0;
        }

         N[E::n_order[0]]     += 0.25 * (1 -     r) * (1 -     s) * (-r - s - 1);
         N[E::n_order[1]]     += 0.25 * (1 +     r) * (1 -     s) * (+r - s - 1);
         N[E::n_order[2]]     += 0.25 * (1 +     r) * (1 +     s) * (+r + s - 1);
         N[E::n_order[3]]     += 0.25 * (1 -     r) * (1 +     s) * (-r + s - 1);
         N[E::n_order[4]]     += 0.5  * (1 - r * r) * (1 -     s);
         N[E::n_order[5]]     += 0.5  * (1 +     r) * (1 - s * s);
         N[E::n_order[6]]     += 0.5  * (1 - r * r) * (1 +     s);
         N[E::n_order[7]]     += 0.5  * (1 -     r) * (1 - s * s);

         dN[E::n_order[0]][0] += 0.25 * (  -     1) * (1 -     s) * (-r - s - 1) - 0.25 * (1 -     r) * (1 -     s);
         dN[E::n_order[1]][0] += 0.25 * (  +     1) * (1 -     s) * (+r - s - 1) + 0.25 * (1 +     r) * (1 -     s);
         dN[E::n_order[2]][0] += 0.25 * (  +     1) * (1 +     s) * (+r + s - 1) + 0.25 * (1 +     r) * (1 +     s);
         dN[E::n_order[3]][0] += 0.25 * (  -     1) * (1 +     s) * (-r + s - 1) - 0.25 * (1 -     r) * (1 +     s);
         dN[E::n_order[4]][0] += 0.5  * (  - 2 * r) * (1 -     s);
         dN[E::n_order[5]][0] += 0.5  * (  +     1) * (1 - s * s);
         dN[E::n_order[6]][0] += 0.5  * (  - 2 * r) * (1 +     s);
         dN[E::n_order[7]][0] += 0.5  * (  -     1) * (1 - s * s);

         dN[E::n_order[0]][1] += 0.25 * (1 -     r) * (  -     1) * (-r - s - 1) - 0.25 * (1 -     r) * (1 -     s);
         dN[E::n_order[1]][1] += 0.25 * (1 +     r) * (  -     1) * (+r - s - 1) - 0.25 * (1 +     r) * (1 -     s);
         dN[E::n_order[2]][1] += 0.25 * (1 +     r) * (  +     1) * (+r + s - 1) + 0.25 * (1 +     r) * (1 +     s);
         dN[E::n_order[3]][1] += 0.25 * (1 -     r) * (  +     1) * (-r + s - 1) + 0.25 * (1 -     r) * (1 +     s);
         dN[E::n_order[4]][1] += 0.5  * (1 - r * r) * (  -     1);
         dN[E::n_order[5]][1] += 0.5  * (1 +     r) * (  - 2 * s);
         dN[E::n_order[6]][1] += 0.5  * (1 - r * r) * (  +     1);
         dN[E::n_order[7]][1] += 0.5  * (1 -     r) * (  - 2 * s);
    }
};

template<typename E> struct GaussPointsDegradedSetter<E, 8, 3>
{
    static void set(double N[], double dN[][3], double r, double s, double t)
    {
        for (size_t n = 0; n < E::nodes; ++n) {
            N[n] = dN[n][0] = dN[n][1] = dN[n][2] = 0;
        }
         N[E::n_order[0]]    += 0.125 * (1 - r) * (1 - s) * (1 - t);
         N[E::n_order[1]]    += 0.125 * (1 + r) * (1 - s) * (1 - t);
         N[E::n_order[2]]    += 0.125 * (1 + r) * (1 + s) * (1 - t);
         N[E::n_order[3]]    += 0.125 * (1 - r) * (1 + s) * (1 - t);
         N[E::n_order[4]]    += 0.125 * (1 - r) * (1 - s) * (1 + t);
         N[E::n_order[5]]    += 0.125 * (1 + r) * (1 - s) * (1 + t);
         N[E::n_order[6]]    += 0.125 * (1 + r) * (1 + s) * (1 + t);
         N[E::n_order[7]]    += 0.125 * (1 - r) * (1 + s) * (1 + t);

        dN[E::n_order[0]][0] += 0.125 * (  - 1) * (1 - s) * (1 - t);
        dN[E::n_order[1]][0] += 0.125 * (  + 1) * (1 - s) * (1 - t);
        dN[E::n_order[2]][0] += 0.125 * (  + 1) * (1 + s) * (1 - t);
        dN[E::n_order[3]][0] += 0.125 * (  - 1) * (1 + s) * (1 - t);
        dN[E::n_order[4]][0] += 0.125 * (  - 1) * (1 - s) * (1 + t);
        dN[E::n_order[5]][0] += 0.125 * (  + 1) * (1 - s) * (1 + t);
        dN[E::n_order[6]][0] += 0.125 * (  + 1) * (1 + s) * (1 + t);
        dN[E::n_order[7]][0] += 0.125 * (  - 1) * (1 + s) * (1 + t);

        dN[E::n_order[0]][1] += 0.125 * (1 - r) * (  - 1) * (1 - t);
        dN[E::n_order[1]][1] += 0.125 * (1 + r) * (  - 1) * (1 - t);
        dN[E::n_order[2]][1] += 0.125 * (1 + r) * (  + 1) * (1 - t);
        dN[E::n_order[3]][1] += 0.125 * (1 - r) * (  + 1) * (1 - t);
        dN[E::n_order[4]][1] += 0.125 * (1 - r) * (  - 1) * (1 + t);
        dN[E::n_order[5]][1] += 0.125 * (1 + r) * (  - 1) * (1 + t);
        dN[E::n_order[6]][1] += 0.125 * (1 + r) * (  + 1) * (1 + t);
        dN[E::n_order[7]][1] += 0.125 * (1 - r) * (  + 1) * (1 + t);

        dN[E::n_order[0]][2] += 0.125 * (1 - r) * (1 - s) * (  - 1);
        dN[E::n_order[1]][2] += 0.125 * (1 + r) * (1 - s) * (  - 1);
        dN[E::n_order[2]][2] += 0.125 * (1 + r) * (1 + s) * (  - 1);
        dN[E::n_order[3]][2] += 0.125 * (1 - r) * (1 + s) * (  - 1);
        dN[E::n_order[4]][2] += 0.125 * (1 - r) * (1 - s) * (  + 1);
        dN[E::n_order[5]][2] += 0.125 * (1 + r) * (1 - s) * (  + 1);
        dN[E::n_order[6]][2] += 0.125 * (1 + r) * (1 + s) * (  + 1);
        dN[E::n_order[7]][2] += 0.125 * (1 - r) * (1 + s) * (  + 1);
    }
};

template<typename E> struct GaussPointsDegradedSetter<E, 20, 3>
{
    static void set(double N[], double dN[][3], double r, double s, double t)
    {
        for (size_t n = 0; n < E::nodes; ++n) {
            N[n] = dN[n][0] = dN[n][1] = dN[n][2] = 0;
        }

         N[E::n_order[ 0]]     += 0.125 * (1 -     r) * (1 -     s) * (1 -     t) * (-r - s - t - 2);
         N[E::n_order[ 1]]     += 0.125 * (1 +     r) * (1 -     s) * (1 -     t) * (+r - s - t - 2);
         N[E::n_order[ 2]]     += 0.125 * (1 +     r) * (1 +     s) * (1 -     t) * (+r + s - t - 2);
         N[E::n_order[ 3]]     += 0.125 * (1 -     r) * (1 +     s) * (1 -     t) * (-r + s - t - 2);
         N[E::n_order[ 4]]     += 0.125 * (1 -     r) * (1 -     s) * (1 +     t) * (-r - s + t - 2);
         N[E::n_order[ 5]]     += 0.125 * (1 +     r) * (1 -     s) * (1 +     t) * (+r - s + t - 2);
         N[E::n_order[ 6]]     += 0.125 * (1 +     r) * (1 +     s) * (1 +     t) * (+r + s + t - 2);
         N[E::n_order[ 7]]     += 0.125 * (1 -     r) * (1 +     s) * (1 +     t) * (-r + s + t - 2);
         N[E::n_order[ 8]]     += 0.25  * (1 - r * r) * (1 -     s) * (1 -     t);
         N[E::n_order[ 9]]     += 0.25  * (1 +     r) * (1 - s * s) * (1 -     t);
         N[E::n_order[10]]     += 0.25  * (1 - r * r) * (1 +     s) * (1 -     t);
         N[E::n_order[11]]     += 0.25  * (1 -     r) * (1 - s * s) * (1 -     t);
         N[E::n_order[12]]     += 0.25  * (1 - r * r) * (1 -     s) * (1 +     t);
         N[E::n_order[13]]     += 0.25  * (1 +     r) * (1 - s * s) * (1 +     t);
         N[E::n_order[14]]     += 0.25  * (1 - r * r) * (1 +     s) * (1 +     t);
         N[E::n_order[15]]     += 0.25  * (1 -     r) * (1 - s * s) * (1 +     t);
         N[E::n_order[16]]     += 0.25  * (1 -     r) * (1 -     s) * (1 - t * t);
         N[E::n_order[17]]     += 0.25  * (1 +     r) * (1 -     s) * (1 - t * t);
         N[E::n_order[18]]     += 0.25  * (1 +     r) * (1 +     s) * (1 - t * t);
         N[E::n_order[19]]     += 0.25  * (1 -     r) * (1 +     s) * (1 - t * t);

         dN[E::n_order[ 0]][0] += 0.125 * (  -     1) * (1 -     s) * (1 -     t) * (-r - s - t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 -     t);
         dN[E::n_order[ 1]][0] += 0.125 * (  +     1) * (1 -     s) * (1 -     t) * (+r - s - t - 2) + 0.125 * (1 +     r) * (1 -     s) * (1 -     t);
         dN[E::n_order[ 2]][0] += 0.125 * (  +     1) * (1 +     s) * (1 -     t) * (+r + s - t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 -     t);
         dN[E::n_order[ 3]][0] += 0.125 * (  -     1) * (1 +     s) * (1 -     t) * (-r + s - t - 2) - 0.125 * (1 -     r) * (1 +     s) * (1 -     t);
         dN[E::n_order[ 4]][0] += 0.125 * (  -     1) * (1 -     s) * (1 +     t) * (-r - s + t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 +     t);
         dN[E::n_order[ 5]][0] += 0.125 * (  +     1) * (1 -     s) * (1 +     t) * (+r - s + t - 2) + 0.125 * (1 +     r) * (1 -     s) * (1 +     t);
         dN[E::n_order[ 6]][0] += 0.125 * (  +     1) * (1 +     s) * (1 +     t) * (+r + s + t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 +     t);
         dN[E::n_order[ 7]][0] += 0.125 * (  -     1) * (1 +     s) * (1 +     t) * (-r + s + t - 2) - 0.125 * (1 -     r) * (1 +     s) * (1 +     t);
         dN[E::n_order[ 8]][0] += 0.25  * (  - 2 * r) * (1 -     s) * (1 -     t);
         dN[E::n_order[ 9]][0] += 0.25  * (  +     1) * (1 - s * s) * (1 -     t);
         dN[E::n_order[10]][0] += 0.25  * (  - 2 * r) * (1 +     s) * (1 -     t);
         dN[E::n_order[11]][0] += 0.25  * (  -     1) * (1 - s * s) * (1 -     t);
         dN[E::n_order[12]][0] += 0.25  * (  - 2 * r) * (1 -     s) * (1 +     t);
         dN[E::n_order[13]][0] += 0.25  * (  +     1) * (1 - s * s) * (1 +     t);
         dN[E::n_order[14]][0] += 0.25  * (  - 2 * r) * (1 +     s) * (1 +     t);
         dN[E::n_order[15]][0] += 0.25  * (  -     1) * (1 - s * s) * (1 +     t);
         dN[E::n_order[16]][0] += 0.25  * (  -     1) * (1 -     s) * (1 - t * t);
         dN[E::n_order[17]][0] += 0.25  * (  +     1) * (1 -     s) * (1 - t * t);
         dN[E::n_order[18]][0] += 0.25  * (  +     1) * (1 +     s) * (1 - t * t);
         dN[E::n_order[19]][0] += 0.25  * (  -     1) * (1 +     s) * (1 - t * t);

         dN[E::n_order[ 0]][1] += 0.125 * (1 -     r) * (  -     1) * (1 -     t) * (-r - s - t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 -     t);
         dN[E::n_order[ 1]][1] += 0.125 * (1 +     r) * (  -     1) * (1 -     t) * (+r - s - t - 2) - 0.125 * (1 +     r) * (1 -     s) * (1 -     t);
         dN[E::n_order[ 2]][1] += 0.125 * (1 +     r) * (  +     1) * (1 -     t) * (+r + s - t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 -     t);
         dN[E::n_order[ 3]][1] += 0.125 * (1 -     r) * (  +     1) * (1 -     t) * (-r + s - t - 2) + 0.125 * (1 -     r) * (1 +     s) * (1 -     t);
         dN[E::n_order[ 4]][1] += 0.125 * (1 -     r) * (  -     1) * (1 +     t) * (-r - s + t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 +     t);
         dN[E::n_order[ 5]][1] += 0.125 * (1 +     r) * (  -     1) * (1 +     t) * (+r - s + t - 2) - 0.125 * (1 +     r) * (1 -     s) * (1 +     t);
         dN[E::n_order[ 6]][1] += 0.125 * (1 +     r) * (  +     1) * (1 +     t) * (+r + s + t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 +     t);
         dN[E::n_order[ 7]][1] += 0.125 * (1 -     r) * (  +     1) * (1 +     t) * (-r + s + t - 2) + 0.125 * (1 -     r) * (1 +     s) * (1 +     t);
         dN[E::n_order[ 8]][1] += 0.25  * (1 - r * r) * (  -     1) * (1 -     t);
         dN[E::n_order[ 9]][1] += 0.25  * (1 +     r) * (  - 2 * s) * (1 -     t);
         dN[E::n_order[10]][1] += 0.25  * (1 - r * r) * (  +     1) * (1 -     t);
         dN[E::n_order[11]][1] += 0.25  * (1 -     r) * (  - 2 * s) * (1 -     t);
         dN[E::n_order[12]][1] += 0.25  * (1 - r * r) * (  -     1) * (1 +     t);
         dN[E::n_order[13]][1] += 0.25  * (1 +     r) * (  - 2 * s) * (1 +     t);
         dN[E::n_order[14]][1] += 0.25  * (1 - r * r) * (  +     1) * (1 +     t);
         dN[E::n_order[15]][1] += 0.25  * (1 -     r) * (  - 2 * s) * (1 +     t);
         dN[E::n_order[16]][1] += 0.25  * (1 -     r) * (  -     1) * (1 - t * t);
         dN[E::n_order[17]][1] += 0.25  * (1 +     r) * (  -     1) * (1 - t * t);
         dN[E::n_order[18]][1] += 0.25  * (1 +     r) * (  +     1) * (1 - t * t);
         dN[E::n_order[19]][1] += 0.25  * (1 -     r) * (  +     1) * (1 - t * t);

         dN[E::n_order[ 0]][2] += 0.125 * (1 -     r) * (1 -     s) * (  -     1) * (-r - s - t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 -     t);
         dN[E::n_order[ 1]][2] += 0.125 * (1 +     r) * (1 -     s) * (  -     1) * (+r - s - t - 2) - 0.125 * (1 +     r) * (1 -     s) * (1 -     t);
         dN[E::n_order[ 2]][2] += 0.125 * (1 +     r) * (1 +     s) * (  -     1) * (+r + s - t - 2) - 0.125 * (1 +     r) * (1 +     s) * (1 -     t);
         dN[E::n_order[ 3]][2] += 0.125 * (1 -     r) * (1 +     s) * (  -     1) * (-r + s - t - 2) - 0.125 * (1 -     r) * (1 +     s) * (1 -     t);
         dN[E::n_order[ 4]][2] += 0.125 * (1 -     r) * (1 -     s) * (  +     1) * (-r - s + t - 2) + 0.125 * (1 -     r) * (1 -     s) * (1 +     t);
         dN[E::n_order[ 5]][2] += 0.125 * (1 +     r) * (1 -     s) * (  +     1) * (+r - s + t - 2) + 0.125 * (1 +     r) * (1 -     s) * (1 +     t);
         dN[E::n_order[ 6]][2] += 0.125 * (1 +     r) * (1 +     s) * (  +     1) * (+r + s + t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 +     t);
         dN[E::n_order[ 7]][2] += 0.125 * (1 -     r) * (1 +     s) * (  +     1) * (-r + s + t - 2) + 0.125 * (1 -     r) * (1 +     s) * (1 +     t);
         dN[E::n_order[ 8]][2] += 0.25  * (1 - r * r) * (1 -     s) * (  -     1);
         dN[E::n_order[ 9]][2] += 0.25  * (1 +     r) * (1 - s * s) * (  -     1);
         dN[E::n_order[10]][2] += 0.25  * (1 - r * r) * (1 +     s) * (  -     1);
         dN[E::n_order[11]][2] += 0.25  * (1 -     r) * (1 - s * s) * (  -     1);
         dN[E::n_order[12]][2] += 0.25  * (1 - r * r) * (1 -     s) * (  +     1);
         dN[E::n_order[13]][2] += 0.25  * (1 +     r) * (1 - s * s) * (  +     1);
         dN[E::n_order[14]][2] += 0.25  * (1 - r * r) * (1 +     s) * (  +     1);
         dN[E::n_order[15]][2] += 0.25  * (1 -     r) * (1 - s * s) * (  +     1);
         dN[E::n_order[16]][2] += 0.25  * (1 -     r) * (1 -     s) * (  - 2 * t);
         dN[E::n_order[17]][2] += 0.25  * (1 +     r) * (1 -     s) * (  - 2 * t);
         dN[E::n_order[18]][2] += 0.25  * (1 +     r) * (1 +     s) * (  - 2 * t);
         dN[E::n_order[19]][2] += 0.25  * (1 -     r) * (1 +     s) * (  - 2 * t);
    }
};

template<typename E, size_t nodes>
struct GaussPointsDegraded<E, nodes, 1>: GaussPointsDegradedSetter<E, nodes, 1> {

    template <typename Element>
    static void simd(Element &element)
    {
        for (size_t r = 0, gp = 0; r < E::order; ++r, ++gp) {
            element.w[gp] = GGaussePoints<E::order>::w[r];
            GaussPointsDegradedSetter<E, nodes, 1>::set(element.N[gp], element.dN[gp], GGaussePoints<E::order>::gps[r]);
        }

        double nn[3] = { -1, 1, 0 };

        for (size_t n = 0; n < nodes; ++n) {
            GaussPointsDegradedSetter<E, nodes, 1>::set(element.NN[n], element.dNN[n], nn[n]);
        }
    }
};

template<typename E, size_t nodes>
struct GaussPointsDegraded<E, nodes, 2>: GaussPointsDegradedSetter<E, nodes, 2> {

    template <typename Element>
    static void simd(Element &element)
    {
        for (size_t r = 0, gp = 0; r < E::order; ++r) {
            for (size_t s = 0; s < E::order; ++s, ++gp) {
                element.w[gp] = GGaussePoints<E::order>::w[r] * GGaussePoints<E::order>::w[s];
                GaussPointsDegradedSetter<E, nodes, 2>::set(element.N[gp], element.dN[gp], GGaussePoints<E::order>::gps[r], GGaussePoints<E::order>::gps[s]);
            }
        }

        double nn[3] = { -1, 1, 0 };

        GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[0]], element.dNN[E::n_order[0]], nn[0], nn[0]);
        GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[1]], element.dNN[E::n_order[1]], nn[1], nn[0]);
        GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[2]], element.dNN[E::n_order[2]], nn[1], nn[1]);
        GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[3]], element.dNN[E::n_order[3]], nn[0], nn[1]);
        if constexpr (nodes == 8) {
            GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[4]], element.dNN[E::n_order[4]], nn[2], nn[0]);
            GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[5]], element.dNN[E::n_order[5]], nn[1], nn[2]);
            GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[6]], element.dNN[E::n_order[6]], nn[2], nn[1]);
            GaussPointsDegradedSetter<E, nodes, 2>::set(element.NN[E::n_order[7]], element.dNN[E::n_order[7]], nn[0], nn[2]);
        }
    }
};

template<typename E, size_t nodes>
struct GaussPointsDegraded<E, nodes, 3>: GaussPointsDegradedSetter<E, nodes, 3> {

    template <typename Element>
    static void simd(Element &element)
    {
        for (size_t r = 0, gp = 0; r < E::order; ++r) {
            for (size_t s = 0; s < E::order; ++s) {
                for (size_t t = 0; t < E::order; ++t, ++gp) {
                    element.w[gp] = GGaussePoints<E::order>::w[r] * GGaussePoints<E::order>::w[s] * GGaussePoints<E::order>::w[t];
                    GaussPointsDegradedSetter<E, nodes, 3>::set(element.N[gp], element.dN[gp], GGaussePoints<E::order>::gps[r], GGaussePoints<E::order>::gps[s], GGaussePoints<E::order>::gps[t]);
                }
            }
        }

        double nn[3] = { -1, 1, 0 };

        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[0]], element.dNN[E::n_order[0]], nn[0], nn[0], nn[0]);
        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[1]], element.dNN[E::n_order[1]], nn[1], nn[0], nn[0]);
        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[2]], element.dNN[E::n_order[2]], nn[1], nn[1], nn[0]);
        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[3]], element.dNN[E::n_order[3]], nn[0], nn[1], nn[0]);
        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[4]], element.dNN[E::n_order[4]], nn[0], nn[0], nn[1]);
        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[5]], element.dNN[E::n_order[5]], nn[1], nn[0], nn[1]);
        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[6]], element.dNN[E::n_order[6]], nn[1], nn[1], nn[1]);
        GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[7]], element.dNN[E::n_order[7]], nn[0], nn[1], nn[1]);

        if constexpr(nodes == 20) {
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[ 8]], element.dNN[E::n_order[ 8]], nn[2], nn[0], nn[0]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[ 9]], element.dNN[E::n_order[ 9]], nn[1], nn[2], nn[0]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[10]], element.dNN[E::n_order[10]], nn[2], nn[1], nn[0]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[11]], element.dNN[E::n_order[11]], nn[0], nn[2], nn[0]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[12]], element.dNN[E::n_order[12]], nn[2], nn[0], nn[1]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[13]], element.dNN[E::n_order[13]], nn[1], nn[2], nn[1]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[14]], element.dNN[E::n_order[14]], nn[2], nn[1], nn[1]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[15]], element.dNN[E::n_order[15]], nn[0], nn[2], nn[1]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[16]], element.dNN[E::n_order[16]], nn[0], nn[0], nn[2]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[17]], element.dNN[E::n_order[17]], nn[1], nn[0], nn[2]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[18]], element.dNN[E::n_order[18]], nn[1], nn[1], nn[2]);
            GaussPointsDegradedSetter<E, nodes, 3>::set(element.NN[E::n_order[19]], element.dNN[E::n_order[19]], nn[0], nn[1], nn[2]);
        }
    }
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_BASEFUNCTIONS_H_ */
