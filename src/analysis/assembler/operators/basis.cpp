
#include "analysis/assembler/module/acoustic.h"
#include "analysis/assembler/module/heattransfer.h"
#include "analysis/assembler/module/structuralmechanics.h"
#include "esinfo/meshinfo.h"
#include "mesh/element.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/boundaryregionstore.h"

namespace espreso {

template<Element::CODE code> void fill(size_t gps, double *N, double *dN, double *w);

template <class Module>
void fill(int code, double *N, double *dN, double *w)
{
	switch ((Element::CODE)code) {
	case Element::CODE::POINT1:    fill<Element::CODE::POINT1>   (Module::NGP::POINT1   , N, dN, w); break;
	case Element::CODE::LINE2:     fill<Element::CODE::LINE2>    (Module::NGP::LINE2    , N, dN, w); break;
	case Element::CODE::TRIANGLE3: fill<Element::CODE::TRIANGLE3>(Module::NGP::TRIANGLE3, N, dN, w); break;
	case Element::CODE::SQUARE4:   fill<Element::CODE::SQUARE4>  (Module::NGP::SQUARE4  , N, dN, w); break;
	case Element::CODE::TETRA4:    fill<Element::CODE::TETRA4>   (Module::NGP::TETRA4   , N, dN, w); break;
	case Element::CODE::PYRAMID5:  fill<Element::CODE::PYRAMID5> (Module::NGP::PYRAMID5 , N, dN, w); break;
	case Element::CODE::PRISMA6:   fill<Element::CODE::PRISMA6>  (Module::NGP::PRISMA6  , N, dN, w); break;
	case Element::CODE::HEXA8:     fill<Element::CODE::HEXA8>    (Module::NGP::HEXA8    , N, dN, w); break;
	case Element::CODE::LINE3:     fill<Element::CODE::LINE3>    (Module::NGP::LINE3    , N, dN, w); break;
	case Element::CODE::TRIANGLE6: fill<Element::CODE::TRIANGLE6>(Module::NGP::TRIANGLE6, N, dN, w); break;
	case Element::CODE::SQUARE8:   fill<Element::CODE::SQUARE8>  (Module::NGP::SQUARE8  , N, dN, w); break;
	case Element::CODE::TETRA10:   fill<Element::CODE::TETRA10>  (Module::NGP::TETRA10  , N, dN, w); break;
	case Element::CODE::PYRAMID13: fill<Element::CODE::PYRAMID13>(Module::NGP::PYRAMID13, N, dN, w); break;
	case Element::CODE::PRISMA15:  fill<Element::CODE::PRISMA15> (Module::NGP::PRISMA15 , N, dN, w); break;
	case Element::CODE::HEXA20:    fill<Element::CODE::HEXA20>   (Module::NGP::HEXA20   , N, dN, w); break;
	default: break;
	}
}

template <class Module>
void _baseFunction(Module &module)
{
	module.integration.N.resize();
	module.integration.dN.resize();
	module.integration.weight.resize();
	{
		int index = 0;
		for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei, ++index) {
			module.integration.N.update[index] = module.integration.dN.update[index] = module.integration.weight.update[index] = 0;

			double *n = (module.integration.N.data->begin() + index)->data();
			double *dn = (module.integration.dN.data->begin() + index)->data();
			double *w = (module.integration.weight.data->begin() + index)->data();

			fill<Module>(ei->code, n, dn, w);

			if (module.settings.simd) {
				esint nodes = Mesh::edata[ei->code].nodes;
				esint gps = Mesh::edata[ei->code].gps;
				esint dim = Mesh::edata[ei->code].dimension;

				for (esint node = nodes - 1; 0 <= node; --node) {
					for (esint gp = gps - 1; 0 <= gp; --gp) {
						for (size_t s = 0; s < SIMD::size; ++s) {
							n[SIMD::size * (gps * node + gp) + s] = n[gps * node + gp];
						}
						for (esint d = dim - 1; 0 <= d; --d) {
							for (size_t s = 0; s < SIMD::size; ++s) {
								dn[SIMD::size * (2 * gps * node + gp * dim + d) + s] = dn[2 * gps * node + gp * dim + d];
							}
						}
					}
				}
				for (esint gp = gps - 1; 0 <= gp; --gp) {
					for (size_t s = 0; s < SIMD::size; ++s) {
						w[SIMD::size * gp + s] = w[gp];
					}
				}
			}
		}
	}

	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension && info::mesh->boundaryRegions[r]->eintervals.size()) {
			module.integration.boundary.N.regions[r].resize();
			module.integration.boundary.dN.regions[r].resize();
			module.integration.boundary.weight.regions[r].resize();

			int index = 0;
			for (auto ei = info::mesh->boundaryRegions[r]->eintervals.begin(); ei != info::mesh->boundaryRegions[r]->eintervals.end(); ++ei, ++index) {
				module.integration.boundary.N.regions[r].update[index] = module.integration.boundary.dN.regions[r].update[index] = module.integration.boundary.weight.regions[r].update[index] = 0;

				double *n = (module.integration.boundary.N.regions[r].data->begin() + index)->data();
				double *dn = (module.integration.boundary.dN.regions[r].data->begin() + index)->data();
				double *w = (module.integration.boundary.weight.regions[r].data->begin() + index)->data();

				fill<Module>(ei->code, n, dn, w);
			}
		}
	}
}

void baseFunction(HeatTransfer &module)
{
	_baseFunction(module);
}

void baseFunction(Acoustic &module)
{
	_baseFunction(module);
}

void baseFunction(StructuralMechanics &module)
{
	_baseFunction(module);
}

template<> void fill<Element::CODE::POINT1>(size_t gps, double *N, double *dN, double *w)
{

}

template<> void fill<Element::CODE::LINE2>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 2;

	switch (gps) {
	case 2: {
		double s[2] = { 1 / sqrt(3), -1 / sqrt(3) };

		for (size_t gp = 0; gp < gps; gp++) {
			w[gp] = 1;

			N[gp * nodes + 0] = (1 - s[gp]) * 0.5;
			N[gp * nodes + 1] = (1 + s[gp]) * 0.5;

			dN[gp * nodes + 0] = -0.5;
			dN[gp * nodes + 1] =  0.5;
		}
	} break;
	}
}

template<> void fill<Element::CODE::TRIANGLE3>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 3;

	switch (gps) {
	case 6: {
		double s[6] = { 0.445948490915965, 0.445948490915965, 0.108103018168070, 0.091576213509771, 0.091576213509771, 0.816847572980459 };
		double t[6] = { 0.445948490915965, 0.108103018168070, 0.445948490915965, 0.091576213509771, 0.816847572980459, 0.091576213509771 };

		w[0] = 0.111690794839005;
		w[1] = 0.111690794839005;
		w[2] = 0.111690794839005;
		w[3] = 0.054975871827661;
		w[4] = 0.054975871827661;
		w[5] = 0.054975871827661;
		for (size_t gp = 0; gp < gps; gp++) {
			N[gp * nodes + 0] = 1 - s[gp] - t[gp];
			N[gp * nodes + 1] = s[gp];
			N[gp * nodes + 2] = t[gp];

			dN[2 * gp * nodes + 0 * nodes + 0] = -1;
			dN[2 * gp * nodes + 0 * nodes + 1] =  1;
			dN[2 * gp * nodes + 0 * nodes + 2] =  0;

			dN[2 * gp * nodes + 1 * nodes + 0] = -1;
			dN[2 * gp * nodes + 1 * nodes + 1] =  0;
			dN[2 * gp * nodes + 1 * nodes + 2] =  1;
		}
	} break;
	}
}

template<> void fill<Element::CODE::SQUARE4>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 4;

	switch (gps) {
	case 4: {
		double CsQ_scale = 0.577350269189626;
		double s[4] = { -CsQ_scale,  CsQ_scale,  CsQ_scale, -CsQ_scale };
		double t[4] = { -CsQ_scale, -CsQ_scale,  CsQ_scale,  CsQ_scale };

		for (size_t gp = 0; gp < gps; gp++) {
			w[gp] = 1;

			N[gp * nodes + 0] = 0.25 * (1 - s[gp]) * (1 - t[gp]);
			N[gp * nodes + 1] = 0.25 * (s[gp] + 1) * (1 - t[gp]);
			N[gp * nodes + 2] = 0.25 * (s[gp] + 1) * (t[gp] + 1);
			N[gp * nodes + 3] = 0.25 * (1 - s[gp]) * (t[gp] + 1);

			dN[2 * gp * nodes + 0 * nodes + 0] = 0.25 * ( t[gp] - 1);
			dN[2 * gp * nodes + 0 * nodes + 1] = 0.25 * (-t[gp] + 1);
			dN[2 * gp * nodes + 0 * nodes + 2] = 0.25 * ( t[gp] + 1);
			dN[2 * gp * nodes + 0 * nodes + 3] = 0.25 * (-t[gp] - 1);

			dN[2 * gp * nodes + 1 * nodes + 0] = 0.25 * ( s[gp] - 1);
			dN[2 * gp * nodes + 1 * nodes + 1] = 0.25 * (-s[gp] - 1);
			dN[2 * gp * nodes + 1 * nodes + 2] = 0.25 * ( s[gp] + 1);
			dN[2 * gp * nodes + 1 * nodes + 3] = 0.25 * (-s[gp] + 1);
		}
	} break;
	}
}

template<> void fill<Element::CODE::TETRA4>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 4;

	switch (gps) {
	case 4: {
		double r[4] = { 0.5854101966249685, 0.1381966011250105, 0.1381966011250105, 0.1381966011250105 };
		double s[4] = { 0.1381966011250105, 0.1381966011250105, 0.1381966011250105, 0.5854101966249685 };
		double t[4] = { 0.1381966011250105, 0.1381966011250105, 0.5854101966249685, 0.1381966011250105 };

		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::PYRAMID5>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 5;

	switch (gps) {
	case 8: {
		double v = 0.577350269189625953;
		double r[8] = {  v,  v,  v,  v, -v, -v, -v, -v };
		double s[8] = { -v, -v,  v,  v, -v, -v,  v,  v };
		double t[8] = { -v,  v, -v,  v, -v,  v, -v,  v };

		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::PRISMA6>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 6;

	switch (gps) {
	case 9: {
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
		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::HEXA8>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 8;

	switch (gps) {
	case 8: {
		double CsQ_scale = 1 / std::sqrt(3);

		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::LINE3>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 3;

	switch (gps) {
	case 3: {
		double s[3] = { -sqrt(3 / 5.0), 0, sqrt(3 / 5.0) };

		w[0] = 5/9.0;
		w[1] = 8/9.0;
		w[2] = 5/9.0;
		for (size_t gp = 0; gp < gps; gp++) {
			N[gp * nodes + 0] = 0.5 * (s[gp] - 1) * s[gp];
			N[gp * nodes + 1] = 0.5 * (s[gp] + 1) * s[gp];
			N[gp * nodes + 2] = 1 - s[gp] * s[gp];

			dN[gp * nodes + 0] = s[gp] - 0.5;
			dN[gp * nodes + 1] = s[gp] + 0.5;
			dN[gp * nodes + 2] = -2 * s[gp];;
		}
	} break;
	}
}

template<> void fill<Element::CODE::TRIANGLE6>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 6;

	switch (gps) {
	case 6: {
		double s[6] = { 0.091576213509771, 0.816847572980459, 0.091576213509771, 0.445948490915965, 0.108103018168070, 0.445948490915965 };
		double t[6] = { 0.091576213509771, 0.091576213509771, 0.816847572980459, 0.445948490915965, 0.445948490915965, 0.108103018168070 };

		w[0] = 0.109951743655322 / 2.0;
		w[1] = 0.109951743655322 / 2.0;
		w[2] = 0.109951743655322 / 2.0;
		w[3] = 0.223381589678011 / 2.0;
		w[4] = 0.223381589678011 / 2.0;
		w[5] = 0.223381589678011 / 2.0;
		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::SQUARE8>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 8;

	switch (gps) {
	case 9: {
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
		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::TETRA10>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 10;

	switch (gps) {
	case 15: {
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
		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::PYRAMID13>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 13;

	switch (gps) {
	case 14: {
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
		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::PRISMA15>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 15;

	switch (gps) {
	case 9: {
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
		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

template<> void fill<Element::CODE::HEXA20>(size_t gps, double *N, double *dN, double *w)
{
	size_t nodes = 20;

	switch (gps) {
	case 8: {
		double v = 0.577350269189625953;
		double _r[8] = {  v,  v,  v,  v, -v, -v, -v, -v };
		double _s[8] = { -v, -v,  v,  v, -v, -v,  v,  v };
		double _t[8] = { -v,  v, -v,  v, -v,  v, -v,  v };

		for (size_t gp = 0; gp < gps; gp++) {
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
	} break;
	}
}

}


