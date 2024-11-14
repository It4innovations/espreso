
#ifndef SRC_ANALYSIS_ASSEMBLER_GENERAL_BASEFUNCTIONS_H_
#define SRC_ANALYSIS_ASSEMBLER_GENERAL_BASEFUNCTIONS_H_

#include "element.h"
#include "esinfo/eslog.h"
#include "mesh/element.h"

#include <array>

namespace espreso {

template <Element::CODE code> struct BaseFunctions;

template <size_t nodes, size_t edim> struct BaseFunctionsDegraded;

template <> struct BaseFunctionsDegraded<8, 3> {
    static void set(double N[8], double dN[8][3], double r, double s, double t, const std::array<int, 8> &node)
    {
        for (int n = 0; n < 8; ++n) {
            N[n] = dN[0][0] = dN[0][1] = dN[0][2] = 0;
        }
         N[node[0]]    += 0.125 * (1 - r) * (1 - s) * (1 - t);
         N[node[1]]    += 0.125 * (1 + r) * (1 - s) * (1 - t);
         N[node[2]]    += 0.125 * (1 + r) * (1 + s) * (1 - t);
         N[node[3]]    += 0.125 * (1 - r) * (1 + s) * (1 - t);
         N[node[4]]    += 0.125 * (1 - r) * (1 - s) * (1 + t);
         N[node[5]]    += 0.125 * (1 + r) * (1 - s) * (1 + t);
         N[node[6]]    += 0.125 * (1 + r) * (1 + s) * (1 + t);
         N[node[7]]    += 0.125 * (1 - r) * (1 + s) * (1 + t);

        dN[node[0]][0] += 0.125 * (  - 1) * (1 - s) * (1 - t);
        dN[node[1]][0] += 0.125 * (  + 1) * (1 - s) * (1 - t);
        dN[node[2]][0] += 0.125 * (  + 1) * (1 + s) * (1 - t);
        dN[node[3]][0] += 0.125 * (  - 1) * (1 + s) * (1 - t);
        dN[node[4]][0] += 0.125 * (  - 1) * (1 - s) * (1 + t);
        dN[node[5]][0] += 0.125 * (  + 1) * (1 - s) * (1 + t);
        dN[node[6]][0] += 0.125 * (  + 1) * (1 + s) * (1 + t);
        dN[node[7]][0] += 0.125 * (  - 1) * (1 + s) * (1 + t);

        dN[node[0]][1] += 0.125 * (1 - r) * (  - 1) * (1 - t);
        dN[node[1]][1] += 0.125 * (1 + r) * (  - 1) * (1 - t);
        dN[node[2]][1] += 0.125 * (1 + r) * (  + 1) * (1 - t);
        dN[node[3]][1] += 0.125 * (1 - r) * (  + 1) * (1 - t);
        dN[node[4]][1] += 0.125 * (1 - r) * (  - 1) * (1 + t);
        dN[node[5]][1] += 0.125 * (1 + r) * (  - 1) * (1 + t);
        dN[node[6]][1] += 0.125 * (1 + r) * (  + 1) * (1 + t);
        dN[node[7]][1] += 0.125 * (1 - r) * (  + 1) * (1 + t);

        dN[node[0]][2] += 0.125 * (1 - r) * (1 - s) * (  - 1);
        dN[node[1]][2] += 0.125 * (1 + r) * (1 - s) * (  - 1);
        dN[node[2]][2] += 0.125 * (1 + r) * (1 + s) * (  - 1);
        dN[node[3]][2] += 0.125 * (1 - r) * (1 + s) * (  - 1);
        dN[node[4]][2] += 0.125 * (1 - r) * (1 - s) * (  + 1);
        dN[node[5]][2] += 0.125 * (1 + r) * (1 - s) * (  + 1);
        dN[node[6]][2] += 0.125 * (1 + r) * (1 + s) * (  + 1);
        dN[node[7]][2] += 0.125 * (1 - r) * (1 + s) * (  + 1);
    }
};

template <> struct BaseFunctionsDegraded<20, 3> {
    static void set(double N[20], double dN[20][3], double r, double s, double t, const std::array<int, 20> &node)
    {
        for (int n = 0; n < 20; ++n) {
            N[n] = dN[0][0] = dN[0][1] = dN[0][2] = 0;
        }

         N[node[ 0]]     += 0.125 * (1 -     r) * (1 -     s) * (1 -     t) * (-r - s - t - 2);
         N[node[ 1]]     += 0.125 * (1 +     r) * (1 -     s) * (1 -     t) * (+r - s - t - 2);
         N[node[ 2]]     += 0.125 * (1 +     r) * (1 +     s) * (1 -     t) * (+r + s - t - 2);
         N[node[ 3]]     += 0.125 * (1 -     r) * (1 +     s) * (1 -     t) * (-r + s - t - 2);
         N[node[ 4]]     += 0.125 * (1 -     r) * (1 -     s) * (1 +     t) * (-r - s + t - 2);
         N[node[ 5]]     += 0.125 * (1 +     r) * (1 -     s) * (1 +     t) * (+r - s + t - 2);
         N[node[ 6]]     += 0.125 * (1 +     r) * (1 +     s) * (1 +     t) * (+r + s + t - 2);
         N[node[ 7]]     += 0.125 * (1 -     r) * (1 +     s) * (1 +     t) * (-r + s + t - 2);
         N[node[ 8]]     += 0.25  * (1 - r * r) * (1 -     s) * (1 -     t);
         N[node[ 9]]     += 0.25  * (1 +     r) * (1 - s * s) * (1 -     t);
         N[node[10]]     += 0.25  * (1 - r * r) * (1 +     s) * (1 -     t);
         N[node[11]]     += 0.25  * (1 -     r) * (1 - s * s) * (1 -     t);
         N[node[12]]     += 0.25  * (1 - r * r) * (1 -     s) * (1 +     t);
         N[node[13]]     += 0.25  * (1 +     r) * (1 - s * s) * (1 +     t);
         N[node[14]]     += 0.25  * (1 - r * r) * (1 +     s) * (1 +     t);
         N[node[15]]     += 0.25  * (1 -     r) * (1 - s * s) * (1 +     t);
         N[node[16]]     += 0.25  * (1 -     r) * (1 -     s) * (1 - t * t);
         N[node[17]]     += 0.25  * (1 +     r) * (1 -     s) * (1 - t * t);
         N[node[18]]     += 0.25  * (1 +     r) * (1 +     s) * (1 - t * t);
         N[node[19]]     += 0.25  * (1 -     r) * (1 +     s) * (1 - t * t);

         dN[node[ 0]][0] += 0.125 * (  -     1) * (1 -     s) * (1 -     t) * (-r - s - t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 -     t);
         dN[node[ 1]][0] += 0.125 * (  +     1) * (1 -     s) * (1 -     t) * (+r - s - t - 2) + 0.125 * (1 +     r) * (1 -     s) * (1 -     t);
         dN[node[ 2]][0] += 0.125 * (  +     1) * (1 +     s) * (1 -     t) * (+r + s - t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 -     t);
         dN[node[ 3]][0] += 0.125 * (  -     1) * (1 +     s) * (1 -     t) * (-r + s - t - 2) - 0.125 * (1 -     r) * (1 +     s) * (1 -     t);
         dN[node[ 4]][0] += 0.125 * (  -     1) * (1 -     s) * (1 +     t) * (-r - s + t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 +     t);
         dN[node[ 5]][0] += 0.125 * (  +     1) * (1 -     s) * (1 +     t) * (+r - s + t - 2) + 0.125 * (1 +     r) * (1 -     s) * (1 +     t);
         dN[node[ 6]][0] += 0.125 * (  +     1) * (1 +     s) * (1 +     t) * (+r + s + t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 +     t);
         dN[node[ 7]][0] += 0.125 * (  -     1) * (1 +     s) * (1 +     t) * (-r + s + t - 2) - 0.125 * (1 -     r) * (1 +     s) * (1 +     t);
         dN[node[ 8]][0] += 0.25  * (  - 2 * r) * (1 -     s) * (1 -     t);
         dN[node[ 9]][0] += 0.25  * (  +     1) * (1 - s * s) * (1 -     t);
         dN[node[10]][0] += 0.25  * (  - 2 * r) * (1 +     s) * (1 -     t);
         dN[node[11]][0] += 0.25  * (  -     1) * (1 - s * s) * (1 -     t);
         dN[node[12]][0] += 0.25  * (  - 2 * r) * (1 -     s) * (1 +     t);
         dN[node[13]][0] += 0.25  * (  +     1) * (1 - s * s) * (1 +     t);
         dN[node[14]][0] += 0.25  * (  - 2 * r) * (1 +     s) * (1 +     t);
         dN[node[15]][0] += 0.25  * (  -     1) * (1 - s * s) * (1 +     t);
         dN[node[16]][0] += 0.25  * (  -     1) * (1 -     s) * (1 - t * t);
         dN[node[17]][0] += 0.25  * (  +     1) * (1 -     s) * (1 - t * t);
         dN[node[18]][0] += 0.25  * (  +     1) * (1 +     s) * (1 - t * t);
         dN[node[19]][0] += 0.25  * (  -     1) * (1 +     s) * (1 - t * t);

         dN[node[ 0]][1] += 0.125 * (1 -     r) * (  -     1) * (1 -     t) * (-r - s - t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 -     t);
         dN[node[ 1]][1] += 0.125 * (1 +     r) * (  -     1) * (1 -     t) * (+r - s - t - 2) - 0.125 * (1 +     r) * (1 -     s) * (1 -     t);
         dN[node[ 2]][1] += 0.125 * (1 +     r) * (  +     1) * (1 -     t) * (+r + s - t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 -     t);
         dN[node[ 3]][1] += 0.125 * (1 -     r) * (  +     1) * (1 -     t) * (-r + s - t - 2) + 0.125 * (1 -     r) * (1 +     s) * (1 -     t);
         dN[node[ 4]][1] += 0.125 * (1 -     r) * (  -     1) * (1 +     t) * (-r - s + t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 +     t);
         dN[node[ 5]][1] += 0.125 * (1 +     r) * (  -     1) * (1 +     t) * (+r - s + t - 2) - 0.125 * (1 +     r) * (1 -     s) * (1 +     t);
         dN[node[ 6]][1] += 0.125 * (1 +     r) * (  +     1) * (1 +     t) * (+r + s + t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 +     t);
         dN[node[ 7]][1] += 0.125 * (1 -     r) * (  +     1) * (1 +     t) * (-r + s + t - 2) + 0.125 * (1 -     r) * (1 +     s) * (1 +     t);
         dN[node[ 8]][1] += 0.25  * (1 - r * r) * (  -     1) * (1 -     t);
         dN[node[ 9]][1] += 0.25  * (1 +     r) * (  - 2 * s) * (1 -     t);
         dN[node[10]][1] += 0.25  * (1 - r * r) * (  +     1) * (1 -     t);
         dN[node[11]][1] += 0.25  * (1 -     r) * (  - 2 * s) * (1 -     t);
         dN[node[12]][1] += 0.25  * (1 - r * r) * (  -     1) * (1 +     t);
         dN[node[13]][1] += 0.25  * (1 +     r) * (  - 2 * s) * (1 +     t);
         dN[node[14]][1] += 0.25  * (1 - r * r) * (  +     1) * (1 +     t);
         dN[node[15]][1] += 0.25  * (1 -     r) * (  - 2 * s) * (1 +     t);
         dN[node[16]][1] += 0.25  * (1 -     r) * (  -     1) * (1 - t * t);
         dN[node[17]][1] += 0.25  * (1 +     r) * (  -     1) * (1 - t * t);
         dN[node[18]][1] += 0.25  * (1 +     r) * (  +     1) * (1 - t * t);
         dN[node[19]][1] += 0.25  * (1 -     r) * (  +     1) * (1 - t * t);

         dN[node[ 0]][2] += 0.125 * (1 -     r) * (1 -     s) * (  -     1) * (-r - s - t - 2) - 0.125 * (1 -     r) * (1 -     s) * (1 -     t);
         dN[node[ 1]][2] += 0.125 * (1 +     r) * (1 -     s) * (  -     1) * (+r - s - t - 2) - 0.125 * (1 +     r) * (1 -     s) * (1 -     t);
         dN[node[ 2]][2] += 0.125 * (1 +     r) * (1 +     s) * (  -     1) * (+r + s - t - 2) - 0.125 * (1 +     r) * (1 +     s) * (1 -     t);
         dN[node[ 3]][2] += 0.125 * (1 -     r) * (1 +     s) * (  -     1) * (-r + s - t - 2) - 0.125 * (1 -     r) * (1 +     s) * (1 -     t);
         dN[node[ 4]][2] += 0.125 * (1 -     r) * (1 -     s) * (  +     1) * (-r - s + t - 2) + 0.125 * (1 -     r) * (1 -     s) * (1 +     t);
         dN[node[ 5]][2] += 0.125 * (1 +     r) * (1 -     s) * (  +     1) * (+r - s + t - 2) + 0.125 * (1 +     r) * (1 -     s) * (1 +     t);
         dN[node[ 6]][2] += 0.125 * (1 +     r) * (1 +     s) * (  +     1) * (+r + s + t - 2) + 0.125 * (1 +     r) * (1 +     s) * (1 +     t);
         dN[node[ 7]][2] += 0.125 * (1 -     r) * (1 +     s) * (  +     1) * (-r + s + t - 2) + 0.125 * (1 -     r) * (1 +     s) * (1 +     t);
         dN[node[ 8]][2] += 0.25  * (1 - r * r) * (1 -     s) * (  -     1);
         dN[node[ 9]][2] += 0.25  * (1 +     r) * (1 - s * s) * (  -     1);
         dN[node[10]][2] += 0.25  * (1 - r * r) * (1 +     s) * (  -     1);
         dN[node[11]][2] += 0.25  * (1 -     r) * (1 - s * s) * (  -     1);
         dN[node[12]][2] += 0.25  * (1 - r * r) * (1 -     s) * (  +     1);
         dN[node[13]][2] += 0.25  * (1 +     r) * (1 - s * s) * (  +     1);
         dN[node[14]][2] += 0.25  * (1 - r * r) * (1 +     s) * (  +     1);
         dN[node[15]][2] += 0.25  * (1 -     r) * (1 - s * s) * (  +     1);
         dN[node[16]][2] += 0.25  * (1 -     r) * (1 -     s) * (  - 2 * t);
         dN[node[17]][2] += 0.25  * (1 +     r) * (1 -     s) * (  - 2 * t);
         dN[node[18]][2] += 0.25  * (1 +     r) * (1 +     s) * (  - 2 * t);
         dN[node[19]][2] += 0.25  * (1 -     r) * (1 +     s) * (  - 2 * t);
    }
};

template <Element::CODE code, size_t gps>
struct GaussPointsX
{
    template <typename Element>
    static void simd(Element &element)
    {
        eslog::error("not-implemented base functions.\n");
    }
};

template<>
struct GaussPointsX<Element::CODE::TETRA4, 4> {
    constexpr static int nodes = 4, edim = 3;

    static void set(double N[nodes], double dN[nodes][edim], double r, double s, double t)
    {
        N[0] = r;
        N[1] = s;
        N[2] = t;
        N[3] = 1.0 - r - s - t;

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
        double r[4] = { 0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105 };
        double s[4] = { 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685 };
        double t[4] = { 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105 };

        for (int gp = 0; gp < 4; ++gp) {
            element.w[gp] = 1.0 / 24.0;
            set(element.N[gp], element.dN[gp], s[gp], t[gp]);
        }
    }
};

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

    constexpr  static std::array<double, 4> gps{ -qa, -qb, qb, qa };
    constexpr  static std::array<double, 4> w  {  wa,  wb, wb, wa };
};

template <> struct GGaussePoints<5> {
    constexpr static double wa = (322. - 13. * sqrt(70.)) / 900.;
    constexpr static double wb = (322. + 13. * sqrt(70.)) / 900.;
    constexpr static double wc = 128 / 225.;
    constexpr static double qa = sqrt(5. + 2. * sqrt(10 / 7.)) / 3.;
    constexpr static double qb = sqrt(5. - 2. * sqrt(10 / 7.)) / 3.;
    constexpr static double qc = 0;

    constexpr  static std::array<double, 5> gps{ -qa, -qb, qc, qb, qa };
    constexpr  static std::array<double, 5> w  {  wa,  wb, wc, wb, wa };
};

template <> struct GGaussePoints<6> {
    constexpr static double wa = 0.171324492379170;
    constexpr static double wb = 0.360761573048139;
    constexpr static double wc = 0.467913934572691;
    constexpr static double qa = 0.932469514203152;
    constexpr static double qb = 0.661209386466265;
    constexpr static double qc = 0.238619186083197;

    constexpr  static std::array<double, 6> gps{ -qa, -qb, -qc, qc, qb, qa };
    constexpr  static std::array<double, 6> w  {  wa,  wb,  wc, wc, wb, wa };
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

    constexpr  static std::array<double, 7> gps{ -qa, -qb, -qc, qd, qc, qb, qa };
    constexpr  static std::array<double, 7> w  {  wa,  wb,  wc, wd, wc, wb, wa };
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

    constexpr  static std::array<double, 8> gps{ -qa, -qb, -qc, -qd, qd, qc, qb, qa };
    constexpr  static std::array<double, 8> w  {  wa,  wb,  wc,  wd, wd, wc, wb, wa };
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
        }; break;
    case 1:
        return gps;
    }
    return -1; // error
}

template<typename E, size_t nodes, size_t edim> struct GaussPointsDegradedSetter;

template<typename E> struct GaussPointsDegradedSetter<E, 2, 1>
{
    static void set(double N[E::nodes], double dN[E::nodes][E::edim], double r)
    {
        for (int n = 0; n < 2; ++n) {
            N[n] = dN[0][0] = 0;
        }
         N[E::n_order[0]]    += 0.5 * (1 - r);
         N[E::n_order[1]]    += 0.5 * (1 + r);

        dN[E::n_order[0]][0] += 0.5 * (  - 1);
        dN[E::n_order[1]][0] += 0.5 * (  + 1);
    }
};

template<typename E> struct GaussPointsDegradedSetter<E, 3, 1>
{
    static void set(double N[E::nodes], double dN[E::nodes][E::edim], double r)
    {
        for (int n = 0; n < 3; ++n) {
            N[n] = dN[0][0] = 0;
        }
         N[E::n_order[0]]    += 0.5 * (1 - r);
         N[E::n_order[1]]    += 0.5 * (1 + r);
         N[E::n_order[2]]    += (1 - r * r);

        dN[E::n_order[0]][0] += 0.5 * (  - 1);
        dN[E::n_order[1]][0] += 0.5 * (  + 1);
        dN[E::n_order[2]][0] += (  - 2 * r);
    }
};

template<typename E> struct GaussPointsDegradedSetter<E, 4, 2>
{
    static void set(double N[E::nodes], double dN[E::nodes][E::edim], double r, double s)
    {
        for (int n = 0; n < 4; ++n) {
            N[n] = dN[0][0] = dN[0][1] = 0;
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
    static void set(double N[E::nodes], double dN[E::nodes][E::edim], double r, double s)
    {
        for (int n = 0; n < 8; ++n) {
            N[n] = dN[0][0] = dN[0][1] = 0;
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
    static void set(double N[E::nodes], double dN[E::nodes][E::edim], double r, double s, double t)
    {
        for (int n = 0; n < 8; ++n) {
            N[n] = dN[0][0] = dN[0][1] = dN[0][2] = 0;
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
    static void set(double N[E::nodes], double dN[E::nodes][E::edim], double r, double s, double t)
    {
        for (int n = 0; n < 20; ++n) {
            N[n] = dN[0][0] = dN[0][1] = dN[0][2] = 0;
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

template<typename E, size_t nodes, size_t edim>
struct GaussPointsDegraded: GaussPointsDegradedSetter<E, nodes, edim> {

    template <typename Element>
    static void simd(Element &element)
    {
        for (size_t r = 0, gp = 0; r < E::order; ++r) {
            for (size_t s = 0; s < E::order; ++s) {
                for (size_t t = 0; t < E::order; ++t, ++gp) {
                    element.w[gp] = GGaussePoints<E::order>::w[r] * GGaussePoints<E::order>::w[s] * GGaussePoints<E::order>::w[t];
                    set(element.N[gp], element.dN[gp], GGaussePoints<E::order>::gps[r] * GGaussePoints<E::order>::gps[s] * GGaussePoints<E::order>::gps[t]);
                }
            }
        }
    }
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::TETRA4, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::TETRA4, ngps>, 8, 3> {
    constexpr static int nodes = 4, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 2, 3, 3, 3, 3 };
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::PYRAMID5, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::PYRAMID5, ngps>, 8, 3> {
    constexpr static int nodes = 5, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 3, 4, 4, 4, 4 };
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::PRISMA6, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::PRISMA6, ngps>, 8, 3> {
    constexpr static int nodes = 6, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 2, 3, 4, 5, 5 };
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::HEXA8, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::HEXA8, ngps>, 8, 3> {
    constexpr static int nodes = 8, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 8> n_order = { 0, 1, 2, 3, 4, 5, 6, 7 };
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::TETRA10, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::TETRA10, ngps>, 20, 3> {
    constexpr static int nodes = 10, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 2, 3, 3, 3, 3, 4, 5, 2, 6, 3, 3, 3, 3, 7, 8, 9, 9 };
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::PYRAMID13, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::PYRAMID13, ngps>, 20, 3> {
    constexpr static int nodes = 13, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, 4, 4, 4, 4, 9, 10, 11, 12 };
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::PRISMA15, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::PRISMA6, ngps>, 20, 3> {
    constexpr static int nodes = 15, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 2, 3, 4, 5, 5, 6, 7, 2, 8, 9, 10, 5, 11, 12, 13, 14, 14 };
};

template <size_t ngps>
struct GaussPointsX<Element::CODE::HEXA20, ngps>: GaussPointsDegraded<GaussPointsX<Element::CODE::HEXA8, ngps>, 20, 3> {
    constexpr static int nodes = 20, gps = ngps, edim = 3, order = getorder<gps, edim>();
    constexpr static std::array<int, 20> n_order = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
};

}

#endif /* SRC_ANALYSIS_ASSEMBLER_GENERAL_BASEFUNCTIONS_H_ */
