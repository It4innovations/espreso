
#ifndef SRC_ANALYSIS_ASSEMBLER_OPERATORS_BASIS_H_
#define SRC_ANALYSIS_ASSEMBLER_OPERATORS_BASIS_H_

#include "analysis/assembler/operator.h"
#include "mesh/element.h"

#include <cmath>
#include <vector>

namespace espreso {

template <Element::CODE code, size_t nodes, size_t gps, size_t edim> struct GaussPoints;

template<>
struct GaussPoints<Element::CODE::LINE2, 2, 2, 1> {

	constexpr static int nodes = 2, gps = 2, edim = 1;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double s[2] = { 1 / sqrt(3), -1 / sqrt(3) };
		for (int gp = 0; gp < gps; gp++) {
			w[gp] = 1;

			N[gp * nodes + 0] = (1 - s[gp]) * 0.5;
			N[gp * nodes + 1] = (1 + s[gp]) * 0.5;

			dN[gp * nodes + 0] = -0.5;
			dN[gp * nodes + 1] =  0.5;
		}
	}
};

template<>
struct GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2> {

	constexpr static int nodes = 3, gps = 6, edim = 2;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

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
			N[gp * nodes + 0] = 1 - s[gp] - t[gp];
			N[gp * nodes + 1] = s[gp];
			N[gp * nodes + 2] = t[gp];

			dN[edim * gp * nodes + 0 * nodes + 0] = -1;
			dN[edim * gp * nodes + 0 * nodes + 1] =  1;
			dN[edim * gp * nodes + 0 * nodes + 2] =  0;

			dN[edim * gp * nodes + 1 * nodes + 0] = -1;
			dN[edim * gp * nodes + 1 * nodes + 1] =  0;
			dN[edim * gp * nodes + 1 * nodes + 2] =  1;
		}
	}
};

template<>
struct GaussPoints<Element::CODE::SQUARE4, 4, 4, 2> {

	constexpr static int nodes = 4, gps = 4, edim = 2;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double CsQ_scale = 0.577350269189626;
		double s[4] = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
		double t[4] = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };

		for (int gp = 0; gp < gps; gp++) {
			w[gp] = 1;

			N[gp * nodes + 0] = 0.25 * (1 - s[gp]) * (1 - t[gp]);
			N[gp * nodes + 1] = 0.25 * (s[gp] + 1) * (1 - t[gp]);
			N[gp * nodes + 2] = 0.25 * (s[gp] + 1) * (t[gp] + 1);
			N[gp * nodes + 3] = 0.25 * (1 - s[gp]) * (t[gp] + 1);

			dN[edim * gp * nodes + 0 * nodes + 0] = 0.25 * ( t[gp] - 1);
			dN[edim * gp * nodes + 0 * nodes + 1] = 0.25 * (-t[gp] + 1);
			dN[edim * gp * nodes + 0 * nodes + 2] = 0.25 * ( t[gp] + 1);
			dN[edim * gp * nodes + 0 * nodes + 3] = 0.25 * (-t[gp] - 1);

			dN[edim * gp * nodes + 1 * nodes + 0] = 0.25 * ( s[gp] - 1);
			dN[edim * gp * nodes + 1 * nodes + 1] = 0.25 * (-s[gp] - 1);
			dN[edim * gp * nodes + 1 * nodes + 2] = 0.25 * ( s[gp] + 1);
			dN[edim * gp * nodes + 1 * nodes + 3] = 0.25 * (-s[gp] + 1);
		}
	}
};

template<>
struct GaussPoints<Element::CODE::TETRA4, 4, 4, 3> {

	constexpr static int nodes = 4, gps = 4, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double r[4] = { 0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105 };
		double s[4] = { 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685 };
		double t[4] = { 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105 };

		for (int gp = 0; gp < gps; gp++) {
			w[gp] = 1.0 / 24.0;

			N[gp * nodes + 0] = r[gp];
			N[gp * nodes + 1] = s[gp];
			N[gp * nodes + 2] = t[gp];
			N[gp * nodes + 3] = 1.0 - r[gp] - s[gp] - t[gp];

			dN[3 * gp * nodes + 0 * nodes + 0] =  1.0;
			dN[3 * gp * nodes + 0 * nodes + 1] =  0.0;
			dN[3 * gp * nodes + 0 * nodes + 2] =  0.0;
			dN[3 * gp * nodes + 0 * nodes + 3] = -1.0;

			dN[3 * gp * nodes + 1 * nodes + 0] =  0.0;
			dN[3 * gp * nodes + 1 * nodes + 1] =  0.0;
			dN[3 * gp * nodes + 1 * nodes + 2] =  1.0;
			dN[3 * gp * nodes + 1 * nodes + 3] = -1.0;

			dN[3 * gp * nodes + 2 * nodes + 0] =  0.0;
			dN[3 * gp * nodes + 2 * nodes + 1] =  1.0;
			dN[3 * gp * nodes + 2 * nodes + 2] =  0.0;
			dN[3 * gp * nodes + 2 * nodes + 3] = -1.0;
		}
	}
};

template<>
struct GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3> {

	constexpr static int nodes = 5, gps = 8, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double v = 0.577350269189625953;
		double r[8] = {  v,  v,  v,  v, -v, -v, -v, -v };
		double s[8] = { -v, -v,  v,  v, -v, -v,  v,  v };
		double t[8] = { -v,  v, -v,  v, -v,  v, -v,  v };

		for (int gp = 0; gp < gps; gp++) {
			w[gp] = 1;

			N[gp * nodes + 0] = 0.125 * ((1 - r[gp]) * (1 - s[gp]) * (1 - t[gp]));
			N[gp * nodes + 1] = 0.125 * ((1 + r[gp]) * (1 - s[gp]) * (1 - t[gp]));
			N[gp * nodes + 2] = 0.125 * ((1 + r[gp]) * (1 + s[gp]) * (1 - t[gp]));
			N[gp * nodes + 3] = 0.125 * ((1 - r[gp]) * (1 + s[gp]) * (1 - t[gp]));
			N[gp * nodes + 4] = 0.125 * ( 4 * (1 + t[gp]));

			dN[3 * gp * nodes + 0 * nodes + 0] = 0.125 * (-(1. - s[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 0 * nodes + 1] = 0.125 * ( (1. - s[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 0 * nodes + 2] = 0.125 * ( (1. + s[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 0 * nodes + 3] = 0.125 * (-(1. + s[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 0 * nodes + 4] = 0;

			dN[3 * gp * nodes + 1 * nodes + 0] = 0.125 * (-(1. - r[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 1 * nodes + 1] = 0.125 * (-(1. + r[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 1 * nodes + 2] = 0.125 * ( (1. + r[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 1 * nodes + 3] = 0.125 * ( (1. - r[gp]) * (1. - t[gp]));
			dN[3 * gp * nodes + 1 * nodes + 4] = 0;

			dN[3 * gp * nodes + 2 * nodes + 0] = 0.125 * (-(1. - r[gp]) * (1. - s[gp]));
			dN[3 * gp * nodes + 2 * nodes + 1] = 0.125 * (-(1. + r[gp]) * (1. - s[gp]));
			dN[3 * gp * nodes + 2 * nodes + 2] = 0.125 * (-(1. + r[gp]) * (1. + s[gp]));
			dN[3 * gp * nodes + 2 * nodes + 3] = 0.125 * (-(1. - r[gp]) * (1. + s[gp]));
			dN[3 * gp * nodes + 2 * nodes + 4] = 0.125 * (4.0);
		}
	}
};

template<>
struct GaussPoints<Element::CODE::PRISMA6, 6, 9, 3> {

	constexpr static int nodes = 6, gps = 9, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

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
			N[gp * nodes + 0] = 0.5 * ((1.0 - t[gp]) * (1.0 - r[gp] - s[gp]));
			N[gp * nodes + 1] = 0.5 * ((1.0 - t[gp]) * r[gp]);
			N[gp * nodes + 2] = 0.5 * ((1.0 - t[gp]) * s[gp]);
			N[gp * nodes + 3] = 0.5 * ((1.0 + t[gp]) * (1.0 - r[gp] - s[gp]));
			N[gp * nodes + 4] = 0.5 * ((1.0 + t[gp]) * r[gp]);
			N[gp * nodes + 5] = 0.5 * ((1.0 + t[gp]) * s[gp]);

			dN[3 * gp * nodes + 0 * nodes + 0] =  t[gp] / 2.0 - 1.0 / 2.0;
			dN[3 * gp * nodes + 0 * nodes + 1] = -t[gp] / 2.0 + 1.0 / 2.0;
			dN[3 * gp * nodes + 0 * nodes + 2] =  0.0;
			dN[3 * gp * nodes + 0 * nodes + 3] = -t[gp] / 2.0 - 1.0 / 2.0;
			dN[3 * gp * nodes + 0 * nodes + 4] =  t[gp] / 2.0 + 1.0 / 2.0;
			dN[3 * gp * nodes + 0 * nodes + 5] =  0;

			dN[3 * gp * nodes + 1 * nodes + 0] =  t[gp] / 2.0 - 1.0 / 2.0;
			dN[3 * gp * nodes + 1 * nodes + 1] =  0.0;
			dN[3 * gp * nodes + 1 * nodes + 2] = -t[gp] / 2.0 + 1.0 / 2.0;
			dN[3 * gp * nodes + 1 * nodes + 3] = -t[gp] / 2.0 - 1.0 / 2.0;
			dN[3 * gp * nodes + 1 * nodes + 4] =  0.0;
			dN[3 * gp * nodes + 1 * nodes + 5] =  t[gp] / 2.0 + 1.0 / 2.0;

			dN[3 * gp * nodes + 2 * nodes + 0] =  r[gp] / 2.0 + s[gp] / 2.0 - 1.0 / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 1] = -r[gp] / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 2] =              - s[gp] / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 3] = -r[gp] / 2.0 - s[gp] / 2.0 + 1.0 / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 4] =  r[gp] / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 5] =                s[gp] / 2.0;
		}
	}
};

template<>
struct GaussPoints<Element::CODE::HEXA8, 8, 8, 3> {

	constexpr static int nodes = 8, gps = 8, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double CsQ_scale = 1 / std::sqrt(3);

		for (int gp = 0; gp < gps; gp++) {
			double r = (gp & 4) ? CsQ_scale : -CsQ_scale;
			double s = (gp & 2) ? CsQ_scale : -CsQ_scale;
			double t = (gp & 1) ? CsQ_scale : -CsQ_scale;

			w[gp] = 1;

			N[gp * nodes + 0] = 0.125 * (1 - r) * (1 - s) * (1 - t);
			N[gp * nodes + 1] = 0.125 * (r + 1) * (1 - s) * (1 - t);
			N[gp * nodes + 2] = 0.125 * (r + 1) * (s + 1) * (1 - t);
			N[gp * nodes + 3] = 0.125 * (1 - r) * (s + 1) * (1 - t);
			N[gp * nodes + 4] = 0.125 * (1 - r) * (1 - s) * (t + 1);
			N[gp * nodes + 5] = 0.125 * (r + 1) * (1 - s) * (t + 1);
			N[gp * nodes + 6] = 0.125 * (r + 1) * (s + 1) * (t + 1);
			N[gp * nodes + 7] = 0.125 * (1 - r) * (s + 1) * (t + 1);

			dN[3 * gp * nodes + 0 * nodes + 0] = 0.125 * (-(1 - s) * (1 - t));
			dN[3 * gp * nodes + 0 * nodes + 1] = 0.125 * ( (1 - s) * (1 - t));
			dN[3 * gp * nodes + 0 * nodes + 2] = 0.125 * ( (1 + s) * (1 - t));
			dN[3 * gp * nodes + 0 * nodes + 3] = 0.125 * (-(1 + s) * (1 - t));
			dN[3 * gp * nodes + 0 * nodes + 4] = 0.125 * (-(1 - s) * (1 + t));
			dN[3 * gp * nodes + 0 * nodes + 5] = 0.125 * ( (1 - s) * (1 + t));
			dN[3 * gp * nodes + 0 * nodes + 6] = 0.125 * ( (1 + s) * (1 + t));
			dN[3 * gp * nodes + 0 * nodes + 7] = 0.125 * (-(1 + s) * (1 + t));

			dN[3 * gp * nodes + 1 * nodes + 0] = 0.125 * (-(1 - r) * (1 - t));
			dN[3 * gp * nodes + 1 * nodes + 1] = 0.125 * (-(1 + r) * (1 - t));
			dN[3 * gp * nodes + 1 * nodes + 2] = 0.125 * ( (1 + r) * (1 - t));
			dN[3 * gp * nodes + 1 * nodes + 3] = 0.125 * ( (1 - r) * (1 - t));
			dN[3 * gp * nodes + 1 * nodes + 4] = 0.125 * (-(1 - r) * (1 + t));
			dN[3 * gp * nodes + 1 * nodes + 5] = 0.125 * (-(1 + r) * (1 + t));
			dN[3 * gp * nodes + 1 * nodes + 6] = 0.125 * ( (1 + r) * (1 + t));
			dN[3 * gp * nodes + 1 * nodes + 7] = 0.125 * ( (1 - r) * (1 + t));

			dN[3 * gp * nodes + 2 * nodes + 0] = 0.125 * (-(1 - r) * (1 - s));
			dN[3 * gp * nodes + 2 * nodes + 1] = 0.125 * (-(1 + r) * (1 - s));
			dN[3 * gp * nodes + 2 * nodes + 2] = 0.125 * (-(1 + r) * (1 + s));
			dN[3 * gp * nodes + 2 * nodes + 3] = 0.125 * (-(1 - r) * (1 + s));
			dN[3 * gp * nodes + 2 * nodes + 4] = 0.125 * ( (1 - r) * (1 - s));
			dN[3 * gp * nodes + 2 * nodes + 5] = 0.125 * ( (1 + r) * (1 - s));
			dN[3 * gp * nodes + 2 * nodes + 6] = 0.125 * ( (1 + r) * (1 + s));
			dN[3 * gp * nodes + 2 * nodes + 7] = 0.125 * ( (1 - r) * (1 + s));
		}
	}
};

template<>
struct GaussPoints<Element::CODE::LINE3, 3, 3, 1> {

	constexpr static int nodes = 3, gps = 3, edim = 1;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double s[3] = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

		w[0] = 5/9.0;
		w[1] = 8/9.0;
		w[2] = 5/9.0;
		for (int gp = 0; gp < gps; gp++) {
			N[gp * nodes + 0] = 0.5 * (s[gp] - 1) * s[gp];
			N[gp * nodes + 1] = 0.5 * (s[gp] + 1) * s[gp];
			N[gp * nodes + 2] = 1 - s[gp] * s[gp];

			dN[gp * nodes + 0] = s[gp] - 0.5;
			dN[gp * nodes + 1] = s[gp] + 0.5;
			dN[gp * nodes + 2] = -2 * s[gp];;
		}
	}
};

template<>
struct GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2> {

	constexpr static int nodes = 6, gps = 6, edim = 2;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

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
			N[gp * nodes + 0] = (1.0 - s[gp] - t[gp]) * (1.0 - 2.0 * (s[gp] + t[gp]));
			N[gp * nodes + 1] = -(s[gp]) * (1.0 - 2.0 * s[gp]);
			N[gp * nodes + 2] = -(t[gp]) * (1.0 - 2.0 * t[gp]);
			N[gp * nodes + 3] = 4.0 * (s[gp]) * (1.0 - s[gp] - t[gp]);
			N[gp * nodes + 4] = 4.0 * (s[gp]) * (t[gp]);
			N[gp * nodes + 5] = 4.0 * (t[gp]) * (1.0 - s[gp] - t[gp]);

			dN[2 * gp * nodes + 0 * nodes + 0] = -3.0 + 4.0 * s[gp] + 4.0 * t[gp];
			dN[2 * gp * nodes + 0 * nodes + 1] = -1.0 + 4.0 * s[gp];
			dN[2 * gp * nodes + 0 * nodes + 2] = 0.0;
			dN[2 * gp * nodes + 0 * nodes + 3] = 4.0 - 8.0 * s[gp] - 4.0 * t[gp];
			dN[2 * gp * nodes + 0 * nodes + 4] = 4.0 * t[gp];
			dN[2 * gp * nodes + 0 * nodes + 5] = -4.0 * t[gp];

			dN[2 * gp * nodes + 1 * nodes + 0] = -3.0 + 4.0 * s[gp] + 4.0 * t[gp];
			dN[2 * gp * nodes + 1 * nodes + 1] = 0.0;
			dN[2 * gp * nodes + 1 * nodes + 2] = -1.0 + 4.0 * t[gp];
			dN[2 * gp * nodes + 1 * nodes + 3] = -4.0 * s[gp];
			dN[2 * gp * nodes + 1 * nodes + 4] = 4.0 * s[gp];
			dN[2 * gp * nodes + 1 * nodes + 5] = 4.0 - 4.0 * s[gp] - 8.0 * t[gp];
		}
	}
};

template<>
struct GaussPoints<Element::CODE::SQUARE8, 8, 9, 2> {

	constexpr static int nodes = 8, gps = 9, edim = 2;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

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
			N[gp * nodes + 0] = -.25 * (s[gp] - 1) * (t[gp] - 1) * (s[gp] + t[gp] + 1);
			N[gp * nodes + 1] =  .25 * (t[gp] - 1) * (-s[gp] * s[gp] + t[gp] * s[gp] + t[gp] + 1);
			N[gp * nodes + 2] =  .25 * (s[gp] + 1) * (t[gp] + 1) * (s[gp] + t[gp] - 1);
			N[gp * nodes + 3] =  .25 * (s[gp] - 1) * (s[gp] - t[gp] + 1) * (t[gp] + 1);
			N[gp * nodes + 4] =  .5  * (s[gp] * s[gp] - 1) * (t[gp] - 1);
			N[gp * nodes + 5] = -.5  * (s[gp] + 1) * (t[gp] * t[gp] - 1);
			N[gp * nodes + 6] = -.5  * (s[gp] * s[gp] - 1) * (t[gp] + 1);
			N[gp * nodes + 7] =  .5  * (s[gp] - 1) * (t[gp] * t[gp] - 1);

			dN[2 * gp * nodes + 0 * nodes + 0] = -((2 * s[gp] + t[gp]) * (t[gp] - 1)) * .25;
			dN[2 * gp * nodes + 0 * nodes + 1] = -((2 * s[gp] - t[gp]) * (t[gp] - 1)) * .25;
			dN[2 * gp * nodes + 0 * nodes + 2] =  ((2 * s[gp] + t[gp]) * (t[gp] + 1)) * .25;
			dN[2 * gp * nodes + 0 * nodes + 3] =  ((2 * s[gp] - t[gp]) * (t[gp] + 1)) * .25;
			dN[2 * gp * nodes + 0 * nodes + 4] = s[gp] * (t[gp] - 1);
			dN[2 * gp * nodes + 0 * nodes + 5] = .5 - t[gp] * t[gp] * .5;
			dN[2 * gp * nodes + 0 * nodes + 6] = -s[gp] * (t[gp] + 1);
			dN[2 * gp * nodes + 0 * nodes + 7] = t[gp] * t[gp] * .5 - .5;

			dN[2 * gp * nodes + 1 * nodes + 0] = -((s[gp] + 2 * t[gp]) * (s[gp] - 1)) * .25;
			dN[2 * gp * nodes + 1 * nodes + 1] = -((s[gp] - 2 * t[gp]) * (s[gp] + 1)) * .25;
			dN[2 * gp * nodes + 1 * nodes + 2] =  ((s[gp] + 2 * t[gp]) * (s[gp] + 1)) * .25;
			dN[2 * gp * nodes + 1 * nodes + 3] =  ((s[gp] - 2 * t[gp]) * (s[gp] - 1)) * .25;
			dN[2 * gp * nodes + 1 * nodes + 4] = s[gp] * s[gp] * .5 - .5;
			dN[2 * gp * nodes + 1 * nodes + 5] = -t[gp] * (s[gp] + 1);
			dN[2 * gp * nodes + 1 * nodes + 6] = .5 - s[gp] * s[gp] * .5;
			dN[2 * gp * nodes + 1 * nodes + 7] = t[gp] * (s[gp] - 1);
		}
	}
};

template<>
struct GaussPoints<Element::CODE::TETRA10, 10, 15, 3> {

	constexpr static int nodes = 10, gps = 15, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

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
			N[gp * nodes + 0] = r[gp] * (2.0 * r[gp] - 1.0);
			N[gp * nodes + 1] = s[gp] * (2.0 * s[gp] - 1.0);
			N[gp * nodes + 2] = t[gp] * (2.0 * t[gp] - 1.0);
			N[gp * nodes + 3] = 2.0 * r[gp] * r[gp] + 4.0 * r[gp] * s[gp] + 4.0 * r[gp] * t[gp] - 3.0 * r[gp] + 2.0* s[gp] * s[gp] + 4.0 * s[gp] * t[gp] - 3.0 * s[gp] + 2.0 * t[gp] * t[gp] - 3.0 * t[gp] + 1.0;
			N[gp * nodes + 4] = 4.0 * r[gp] * s[gp];
			N[gp * nodes + 5] = 4.0 * s[gp] * t[gp];
			N[gp * nodes + 6] = 4.0 * r[gp] * t[gp];
			N[gp * nodes + 7] = r[gp] * (-4.0 * r[gp] - 4.0 * s[gp] - 4.0 * t[gp] + 4.0);
			N[gp * nodes + 8] = s[gp] * (-4.0 * r[gp] - 4.0 * s[gp] - 4.0 * t[gp] + 4.0);
			N[gp * nodes + 9] = t[gp] * (-4.0 * r[gp] - 4.0 * s[gp] - 4.0 * t[gp] + 4.0);

			dN[3 * gp * nodes + 0 * nodes + 0] = 4.0 * r[gp] - 1.0;
			dN[3 * gp * nodes + 0 * nodes + 1] = 0;
			dN[3 * gp * nodes + 0 * nodes + 2] = 0;
			dN[3 * gp * nodes + 0 * nodes + 3] = 4.0 * r[gp] + 4.0 * s[gp] + 4.0 * t[gp] - 3.0;
			dN[3 * gp * nodes + 0 * nodes + 4] = 4.0 * s[gp];
			dN[3 * gp * nodes + 0 * nodes + 5] = 0;
			dN[3 * gp * nodes + 0 * nodes + 6] = 4.0 * t[gp];
			dN[3 * gp * nodes + 0 * nodes + 7] = -8.0 * r[gp] - 4.0 * s[gp] - 4.0 * t[gp] + 4.0;
			dN[3 * gp * nodes + 0 * nodes + 8] = -4.0 * s[gp];
			dN[3 * gp * nodes + 0 * nodes + 9] = -4.0 * t[gp];

			dN[3 * gp * nodes + 1 * nodes + 0] = 0;
			dN[3 * gp * nodes + 1 * nodes + 1] = 0;
			dN[3 * gp * nodes + 1 * nodes + 2] = 4.0 * t[gp] - 1.0;
			dN[3 * gp * nodes + 1 * nodes + 3] = 4.0 * r[gp] + 4.0 * s[gp] + 4.0* t[gp]  - 3.0;
			dN[3 * gp * nodes + 1 * nodes + 4] = 0;
			dN[3 * gp * nodes + 1 * nodes + 5] = 4.0 * s[gp];
			dN[3 * gp * nodes + 1 * nodes + 6] = 4.0 * r[gp];
			dN[3 * gp * nodes + 1 * nodes + 7] = -4.0 * r[gp];
			dN[3 * gp * nodes + 1 * nodes + 8] = -4.0 * s[gp];
			dN[3 * gp * nodes + 1 * nodes + 9] = -4.0 * r[gp] - 4.0 * s[gp] - 8.0 * t[gp] + 4.0;

			dN[3 * gp * nodes + 2 * nodes + 0] = 0;
			dN[3 * gp * nodes + 2 * nodes + 1] = 4.0 * s[gp] - 1.0;
			dN[3 * gp * nodes + 2 * nodes + 2] = 0 ;
			dN[3 * gp * nodes + 2 * nodes + 3] = 4.0 * r[gp] + 4.0 * s[gp] + 4.0 * t[gp] - 3.0;
			dN[3 * gp * nodes + 2 * nodes + 4] = 4.0 * r[gp];
			dN[3 * gp * nodes + 2 * nodes + 5] = 4.0 * t[gp];
			dN[3 * gp * nodes + 2 * nodes + 6] = 0;
			dN[3 * gp * nodes + 2 * nodes + 7] = -4.0 * r[gp];
			dN[3 * gp * nodes + 2 * nodes + 8] = -4.0 * r[gp] - 8.0 * s[gp] - 4.0 * t[gp] + 4.0;
			dN[3 * gp * nodes + 2 * nodes + 9] = -4.0 * t[gp];
		}
	}
};

template<>
struct GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3> {

	constexpr static int nodes = 13, gps = 14, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double v1 = 0.758786910639329015;
		double v2 = 0.795822425754222018;
		double v3 = 0;
		double _r[14] = { -v1,  v1,  v1, -v1, -v1,  v1,  v1, -v1,  v3,  v3,  v2, v3, -v2, v3 };
		double _s[14] = { -v1, -v1,  v1,  v1, -v1, -v1,  v1,  v1,  v3, -v2,  v3, v2,  v3, v3 };
		double _t[14] = { -v1, -v1, -v1, -v1,  v1,  v1,  v1,  v1, -v2,  v3,  v3, v3,  v3, v2 };

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
			double r = _r[gp];
			double s = _s[gp];
			double t = _t[gp];

			N[gp * nodes +  0] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 - s) * (-1.0 - (0.5 * (1.0 - t)) * r - (0.5 * (1.0 - t)) * s));
			N[gp * nodes +  1] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 - s) * (-1.0 + (0.5 * (1 - t))   * r - (0.5 * (1.0 - t)) * s));
			N[gp * nodes +  2] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 + r) * (1.0 + s) * (-1.0 + (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
			N[gp * nodes +  3] = ((0.5 * (1.0 - t)) / 4.0) * ((1.0 - r) * (1.0 + s) * (-1.0 - (0.5 * (1.0 - t)) * r + (0.5 * (1.0 - t)) * s));
			N[gp * nodes +  4] = (1.0 - (0.5 * (1.0 - t))) * (1.0 - 2.0 * (0.5 * (1.0 - t)));
			N[gp * nodes +  5] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - s) * (1.0 - r * r);
			N[gp * nodes +  6] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + r) * (1.0 - s * s);
			N[gp * nodes +  7] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 + s) * (1.0 - r * r);
			N[gp * nodes +  8] = (((0.5 * (1.0 - t)) * (0.5 * (1.0 - t))) / 2.0) * (1.0 - r) * (1.0 - s * s);
			N[gp * nodes +  9] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r - s + r * s);
			N[gp * nodes + 10] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r - s - r * s);
			N[gp * nodes + 11] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 + r + s + r * s);
			N[gp * nodes + 12] = (0.5 * (1.0 - t)) * (1.0 - (0.5 * (1.0 - t))) * (1.0 - r + s - r * s);

			dN[3 * gp * nodes + 0 * nodes +  0] = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  1] = -(t / 8.0 - 1.0 / 8.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0) - (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  2] =  (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0) + (t / 2.0 - 0.5) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
			dN[3 * gp * nodes + 0 * nodes +  3] =  (t / 2.0 - 0.5) * (s + 1.0) * (r - 1.0) * (t / 8.0 - 1.0 / 8.0) - (t / 8.0 - 1.0 / 8.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0);
			dN[3 * gp * nodes + 0 * nodes +  4] =  0.0;
			dN[3 * gp * nodes + 0 * nodes +  5] =  r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  6] = -((s * s - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes +  7] = -r * ((t / 2.0 - 0.5) * (t / 2.0 - 0.5)) * (s + 1.0);
			dN[3 * gp * nodes + 0 * nodes +  8] =  ((s * s - 1) * (t / 2.0 - 0.5) * (t / 2.0 - 0.5)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes +  9] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
			dN[3 * gp * nodes + 0 * nodes + 10] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s - 1.0);
			dN[3 * gp * nodes + 0 * nodes + 11] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);
			dN[3 * gp * nodes + 0 * nodes + 12] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (s + 1.0);

			dN[3 * gp * nodes + 1 * nodes +  0] = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) - 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
			dN[3 * gp * nodes + 1 * nodes +  1] =  (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0) - (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) - s * (t / 2.0 - 1.0 / 2.0) + 1.0);
			dN[3 * gp * nodes + 1 * nodes +  2] =  (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (r * (t / 2.0 - 1.0 / 2.0) + s * (t / 2.0 - 1.0 / 2.0) + 1.0) + (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
			dN[3 * gp * nodes + 1 * nodes +  3] = -(t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s * (t / 2.0 - 1.0 / 2.0) - r * (t / 2.0 - 1.0 / 2.0) + 1.0) - (t / 2.0 - 1.0 / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0);
			dN[3 * gp * nodes + 1 * nodes +  4] =  0.0;
			dN[3 * gp * nodes + 1 * nodes +  5] =  ((r * r - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes +  6] = -s * (t / 2.0 - 1.0 / 2.0) * (t / 2.0 - 0.5) * (r + 1.0);
			dN[3 * gp * nodes + 1 * nodes +  7] = -((r * r - 1.0) * (t / 2.0 - 0.5) * (t / 2.0 - 1.0 / 2.0)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes +  8] =  s * (t / 2.0 - 0.5) * (t / 2.0 - 0.5) * (r - 1.0);
			dN[3 * gp * nodes + 1 * nodes +  9] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);
			dN[3 * gp * nodes + 1 * nodes + 10] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
			dN[3 * gp * nodes + 1 * nodes + 11] = -(t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r + 1.0);
			dN[3 * gp * nodes + 1 * nodes + 12] =  (t / 2.0 - 0.5) * (t / 2.0 + 0.5) * (r - 1.0);

			dN[3 * gp * nodes + 2 * nodes +  0] = -((r - 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) - 1.0)) / 8.0 - (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s - 1.0);
			dN[3 * gp * nodes + 2 * nodes +  1] = -((r + 1.0) * (s - 1.0) * (r * (t / 2.0 - 0.5) - s * (t / 2.0 - 0.5) + 1.0)) / 8.0 - (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s - 1.0);
			dN[3 * gp * nodes + 2 * nodes +  2] =  ((r + 1.0) * (s + 1.0) * (r * (t / 2.0 - 0.5) + s * (t / 2.0 - 0.5) + 1.0)) / 8.0 + (r / 2.0 + s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r + 1.0) * (s + 1.0);
			dN[3 * gp * nodes + 2 * nodes +  3] =  (r / 2.0 - s / 2.0) * (t / 8.0 - 1.0 / 8.0) * (r - 1.0) * (s + 1.0) - ((r - 1.0) * (s + 1.0) * (s * (t / 2.0 - 0.5) - r * (t / 2.0 - 0.5) + 1.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  4] =  t + 0.5;
			dN[3 * gp * nodes + 2 * nodes +  5] =  ((r * r - 1.0) * (t / 2.0 - 0.5) * (s - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  6] = -((s * s - 1.0) * (t / 2.0 - 0.5) * (r + 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  7] = -((r * r - 1.0) * (t / 2.0 - 0.5) * (s + 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  8] =  ((s * s - 1.0) * (t / 2.0 - 0.5) * (r - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  9] =  ((t / 2.0 - 0.5) * (r + s - r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r + s - r * s - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 10] = -((t / 2.0 - 0.5) * (r - s - r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r - s - r * s + 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 11] = -((t / 2.0 - 0.5) * (r + s + r * s + 1.0)) / 2.0 - ((t / 2.0 + 0.5) * (r + s + r * s + 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 12] =  ((t / 2.0 - 0.5) * (r - s + r * s - 1.0)) / 2.0 + ((t / 2.0 + 0.5) * (r - s + r * s - 1.0)) / 2.0;
		}
	}
};

template<>
struct GaussPoints<Element::CODE::PRISMA15, 15, 9, 3> {

	constexpr static int nodes = 15, gps = 9, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double v1 = 1.0 / 6.0;
		double v2 = 4.0 / 6.0;
		double v3 = sqrt(3.0 / 5.0);
		double v4 = 0.0;
		double _r[9] = {  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2,  v1 };
		double _s[9] = {  v1,  v1,  v2,  v1,  v1,  v2,  v1,  v1,  v2 };
		double _t[9] = { -v3, -v3, -v3,  v4,  v4,  v4,  v3,  v3,  v3 };

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
			double r = _r[gp];
			double s = _s[gp];
			double t = _t[gp];

			N[gp * nodes +  0] = -(1.0 - r - s) * (1.0 - t) * (2.0 * r + 2.0 * s + t) / 2.0;
			N[gp * nodes +  1] = r * (1.0 - t) * (2.0 * r - t - 2.0) / 2.0;
			N[gp * nodes +  2] = s * (1.0 - t) * (2.0 * s - t - 2.0) / 2.0;
			N[gp * nodes +  3] = -(1.0 - r - s) * (1.0 + t) * (2.0 * r + 2.0 * s - t) / 2.0;
			N[gp * nodes +  4] = r * (t + 1.0) * (2.0 * r + t - 2.0) / 2.0;
			N[gp * nodes +  5] = s * (t + 1.0) * (2.0 * s + t - 2.0) / 2.0;
			N[gp * nodes +  6] = 2.0 * r * (1.0 - r - s) * (1.0 - t);
			N[gp * nodes +  7] = 2.0 * r * s * (1.0 - t);
			N[gp * nodes +  8] = 2.0 * s * (1.0 - r - s) * (1.0 - t);
			N[gp * nodes +  9] = 2.0 * r * (1.0 - r - s) * (1.0 + t);
			N[gp * nodes + 10] = 2.0 * r * s * (1.0 + t);
			N[gp * nodes + 11] = 2.0 * s * (1.0 - r - s) * (1.0 + t);
			N[gp * nodes + 12] = (1.0 - r - s) * (1.0 - t * t);
			N[gp * nodes + 13] = r * (1.0 - t * t);
			N[gp * nodes + 14] = s * (1.0 - t * t);

			dN[3 * gp * nodes + 0 * nodes +  0] = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes +  1] = ((t - 1.0) * (t - 2.0 * r + 2.0)) / 2.0 - r * (t - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  2] = 0.0;
			dN[3 * gp * nodes + 0 * nodes +  3] = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  4] = r * (t + 1.0) + ((t + 1.0) * (2.0 * r + t - 2.0)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes +  5] = 0.0;
			dN[3 * gp * nodes + 0 * nodes +  6] = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * r * (t - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  7] = (-2.0) * s * (t - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  8] = 2.0 * s * (t - 1.0);
			dN[3 * gp * nodes + 0 * nodes +  9] = -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * r * (t + 1.0);
			dN[3 * gp * nodes + 0 * nodes + 10] =  2.0 * s * (t + 1.0);
			dN[3 * gp * nodes + 0 * nodes + 11] =  -2.0 * s * (t + 1.0);
			dN[3 * gp * nodes + 0 * nodes + 12] =  t * t - 1.0;
			dN[3 * gp * nodes + 0 * nodes + 13] =  1.0 - t * t;
			dN[3 * gp * nodes + 0 * nodes + 14] =  0.0;

			dN[3 * gp * nodes + 1 * nodes +  0] = -(t - 1.0) * (r + s - 1.0) - ((t - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes +  1] = 0.0;
			dN[3 * gp * nodes + 1 * nodes +  2] = ((t - 1.0) * (t - 2.0 * s + 2.0)) / 2.0 - s * (t - 1.0);
			dN[3 * gp * nodes + 1 * nodes +  3] = ((t + 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 + (t + 1.0) * (r + s - 1.0);
			dN[3 * gp * nodes + 1 * nodes +  4] = 0.0;
			dN[3 * gp * nodes + 1 * nodes +  5] = s * (t + 1.0) + ((t + 1.0) * (2.0 * s + t - 2.0)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes +  6] = 2.0 * r * (t - 1.0);
			dN[3 * gp * nodes + 1 * nodes +  7] = (-2.0) * r * (t - 1.0);
			dN[3 * gp * nodes + 1 * nodes +  8] = 2.0 * (t - 1.0) * (r + s - 1.0) + 2.0 * s * (t - 1.0);
			dN[3 * gp * nodes + 1 * nodes +  9] = (-2.0) * r * (t + 1.0);
			dN[3 * gp * nodes + 1 * nodes + 10] =  2.0 * r * (t + 1.0);
			dN[3 * gp * nodes + 1 * nodes + 11] =  -2.0 * (t + 1.0) * (r + s - 1.0) - 2.0 * s * (t + 1.0);
			dN[3 * gp * nodes + 1 * nodes + 12] =  t * t - 1.0;
			dN[3 * gp * nodes + 1 * nodes + 13] =  0.0;
			dN[3 * gp * nodes + 1 * nodes + 14] =  1.0 - t * t;

			dN[3 * gp * nodes + 2 * nodes +  0] = -((r + s - 1.0) * (2.0 * r + 2.0 * s + t)) / 2.0 - ((t - 1.0) * (r + s - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  1] = (r * (t - 2.0 * r + 2.0)) / 2.0 + (r * (t - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  2] = (s * (t - 2.0 * s + 2.0)) / 2.0 + (s * (t - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  3] = ((r + s - 1.0) * (2.0 * r + 2.0 * s - t)) / 2.0 - ((t + 1.0) * (r + s - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  4] = (r * (2.0 * r + t - 2.0)) / 2.0 + (r * (t + 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  5] = (s * (2.0 * s + t - 2.0)) / 2.0 + (s * (t + 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes +  6] = 2.0 * r * (r + s - 1.0);
			dN[3 * gp * nodes + 2 * nodes +  7] = (-2.0) * r * s;
			dN[3 * gp * nodes + 2 * nodes +  8] = 2.0 * s * (r + s - 1.0);
			dN[3 * gp * nodes + 2 * nodes +  9] = (-2.0) * r * (r + s - 1.0);
			dN[3 * gp * nodes + 2 * nodes + 10] =  2.0 * r * s;
			dN[3 * gp * nodes + 2 * nodes + 11] =  (-2.0) * s * (r + s - 1.0);
			dN[3 * gp * nodes + 2 * nodes + 12] =  2.0 * t * (r + s - 1.0);
			dN[3 * gp * nodes + 2 * nodes + 13] =  (-2.0) * r * t;
			dN[3 * gp * nodes + 2 * nodes + 14] =  (-2.0) * s * t;
		}
	}
};

template<>
struct GaussPoints<Element::CODE::HEXA20, 20, 8, 3> {

	constexpr static int nodes = 20, gps = 8, edim = 3;
	static double w[gps], N[gps * nodes], dN[gps * nodes * edim];

	static void set()
	{
		double v = 0.577350269189625953;
		double _r[8] = {  v,  v,  v,  v, -v, -v, -v, -v };
		double _s[8] = { -v, -v,  v,  v, -v, -v,  v,  v };
		double _t[8] = { -v,  v, -v,  v, -v,  v, -v,  v };

		for (int gp = 0; gp < gps; gp++) {
			double r = _r[gp];
			double s = _s[gp];
			double t = _t[gp];

			w[gp] = 1;

			N[gp * nodes +  0] = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 - t) * (-r - s - t - 2.0));
			N[gp * nodes +  1] = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 - t) * ( r - s - t - 2.0));
			N[gp * nodes +  2] = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 - t) * ( r + s - t - 2.0));
			N[gp * nodes +  3] = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 - t) * (-r + s - t - 2.0));
			N[gp * nodes +  4] = 0.125 * ((1.0 - r) * (1.0 - s) * (1.0 + t) * (-r - s + t - 2.0));
			N[gp * nodes +  5] = 0.125 * ((1.0 + r) * (1.0 - s) * (1.0 + t) * ( r - s + t - 2.0));
			N[gp * nodes +  6] = 0.125 * ((1.0 + r) * (1.0 + s) * (1.0 + t) * ( r + s + t - 2.0));
			N[gp * nodes +  7] = 0.125 * ((1.0 - r) * (1.0 + s) * (1.0 + t) * (-r + s + t - 2.0));
			N[gp * nodes +  8] = 0.25 * ((1.0 - r * r) * (1.0 - s) * (1.0 - t));
			N[gp * nodes +  9] = 0.25 * ((1.0 + r) * (1.0 - s * s) * (1.0 - t));
			N[gp * nodes + 10] = 0.25 * ((1.0 - r * r) * (1.0 + s) * (1.0 - t));
			N[gp * nodes + 11] = 0.25 * ((1.0 - r) * (1.0 - s * s) * (1.0 - t));
			N[gp * nodes + 12] = 0.25 * ((1.0 - r * r) * (1.0 - s) * (1.0 + t));
			N[gp * nodes + 13] = 0.25 * ((1.0 + r) * (1.0 - s * s) * (1.0 + t));
			N[gp * nodes + 14] = 0.25 * ((1.0 - r * r) * (1.0 + s) * (1.0 + t));
			N[gp * nodes + 15] = 0.25 * ((1.0 - r) * (1.0 - s * s) * (1.0 + t));
			N[gp * nodes + 16] = 0.25 * ((1.0 - r) * (1.0 - s) * (1.0 - t * t));
			N[gp * nodes + 17] = 0.25 * ((1.0 + r) * (1.0 - s) * (1.0 - t * t));
			N[gp * nodes + 18] = 0.25 * ((1.0 + r) * (1.0 + s) * (1.0 - t * t));
			N[gp * nodes + 19] = 0.25 * ((1.0 - r) * (1.0 + s) * (1.0 - t * t));

			dN[3 * gp * nodes + 0 * nodes +  0] =  ((s - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  1] =  ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0 - ((s - 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  2] = -((s + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  3] = -((s + 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  4] = -((s - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  5] = -((s - 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  6] =  ((s + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  7] =  ((s + 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 + ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 0 * nodes +  8] =  -(r * (s - 1.0) * (t - 1.0)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes +  9] =   ((s * s - 1.0) * (t - 1.0)) / 4.0;
			dN[3 * gp * nodes + 0 * nodes + 10] =   (r * (s + 1.0) * (t - 1.0)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes + 11] =  -((s * s - 1.0) * (t - 1.0)) / 4.0;
			dN[3 * gp * nodes + 0 * nodes + 12] =   (r * (s - 1.0) * (t + 1.0)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes + 13] =  -((s * s - 1.0) * (t + 1.0)) / 4.0;
			dN[3 * gp * nodes + 0 * nodes + 14] =  -(r * (s + 1.0) * (t + 1.0)) / 2.0;
			dN[3 * gp * nodes + 0 * nodes + 15] =   ((s * s - 1.0) * (t + 1.0)) / 4.0;
			dN[3 * gp * nodes + 0 * nodes + 16] =  -((t * t - 1.0) * (s - 1.0)) / 4.0;
			dN[3 * gp * nodes + 0 * nodes + 17] =   ((t * t - 1.0) * (s - 1.0)) / 4.0;
			dN[3 * gp * nodes + 0 * nodes + 18] =  -((t * t - 1.0) * (s + 1.0)) / 4.0;
			dN[3 * gp * nodes + 0 * nodes + 19] =   ((t * t - 1.0) * (s + 1.0)) / 4.0;

			dN[3 * gp * nodes + 1 * nodes +  0] =  ((r - 1.0) * (t - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  1] = -((r + 1.0) * (t - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  2] = -((r + 1.0) * (t - 1.0) * (r + s - t - 2.0)) / 8.0 - ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  3] =  ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r - 1.0) * (t - 1.0) * (r - s + t + 2.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  4] = -((r - 1.0) * (t + 1.0) * (r + s - t + 2.0)) / 8.0 - ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  5] =  ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r + 1.0) * (t + 1.0) * (r - s + t - 2.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  6] =  ((r + 1.0) * (t + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  7] =  ((r - 1.0) * (t + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 1 * nodes +  8] =  -((r * r - 1.0) * (t - 1.0)) / 4.0;
			dN[3 * gp * nodes + 1 * nodes +  9] =   (s * (r + 1.0) * (t - 1.0)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes + 10] =   ((r * r - 1.0) * (t - 1.0)) / 4.0;
			dN[3 * gp * nodes + 1 * nodes + 11] =  -(s * (r - 1.0) * (t - 1.0)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes + 12] =   ((r * r - 1.0) * (t + 1.0)) / 4.0;
			dN[3 * gp * nodes + 1 * nodes + 13] =  -(s * (r + 1.0) * (t + 1.0)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes + 14] =  -((r * r - 1.0) * (t + 1.0)) / 4.0;
			dN[3 * gp * nodes + 1 * nodes + 15] =   (s * (r - 1.0) * (t + 1.0)) / 2.0;
			dN[3 * gp * nodes + 1 * nodes + 16] =  -((t * t - 1.0) * (r - 1.0)) / 4.0;
			dN[3 * gp * nodes + 1 * nodes + 17] =   ((t * t - 1.0) * (r + 1.0)) / 4.0;
			dN[3 * gp * nodes + 1 * nodes + 18] =  -((t * t - 1.0) * (r + 1.0)) / 4.0;
			dN[3 * gp * nodes + 1 * nodes + 19] =   ((t * t - 1.0) * (r - 1.0)) / 4.0;

			dN[3 * gp * nodes + 2 * nodes +  0] =  ((r - 1.0) * (s - 1.0) * (r + s + t + 2.0)) / 8.0 + ((r - 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  1] = -((r + 1.0) * (s - 1.0) * (s - r + t + 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  2] =  ((r + 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0 - ((r + 1.0) * (s + 1.0) * (r + s - t - 2.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  3] = -((r - 1.0) * (s + 1.0) * (r - s + t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t - 1.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  4] =  ((r - 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0 - ((r - 1.0) * (s - 1.0) * (r + s - t + 2.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  5] = -((r + 1.0) * (s - 1.0) * (r - s + t - 2.0)) / 8.0 - ((r + 1.0) * (s - 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  6] =  ((r + 1.0) * (s + 1.0) * (r + s + t - 2.0)) / 8.0 + ((r + 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  7] =  ((r - 1.0) * (s + 1.0) * (r - s - t + 2.0)) / 8.0 - ((r - 1.0) * (s + 1.0) *         (t + 1.0)) / 8.0;
			dN[3 * gp * nodes + 2 * nodes +  8] =  -((r * r - 1.0) * (s - 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes +  9] =   ((s * s - 1.0) * (r + 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes + 10] =   ((r * r - 1.0) * (s + 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes + 11] =  -((s * s - 1.0) * (r - 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes + 12] =   ((r * r - 1.0) * (s - 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes + 13] =  -((s * s - 1.0) * (r + 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes + 14] =  -((r * r - 1.0) * (s + 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes + 15] =   ((s * s - 1.0) * (r - 1.0)) / 4.0;
			dN[3 * gp * nodes + 2 * nodes + 16] =  -(t * (r - 1.0) * (s - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 17] =   (t * (r + 1.0) * (s - 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 18] =  -(t * (r + 1.0) * (s + 1.0)) / 2.0;
			dN[3 * gp * nodes + 2 * nodes + 19] =   (t * (r - 1.0) * (s + 1.0)) / 2.0;
		}
	}

};

template <Element::CODE code, size_t nodes, size_t gps, size_t edim, size_t etype, class Physics>
struct Basis: ActionOperator, Physics {
	const char* name() const { return "Basis"; }

	Basis()
	{
		action = Action::ASSEMBLE | Action::REASSEMBLE | Action::SOLUTION;
	}

	void simd(typename Physics::Element &element)
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
};


}

#endif /* SRC_ANALYSIS_ASSEMBLER_OPERATORS_BASIS_H_ */
