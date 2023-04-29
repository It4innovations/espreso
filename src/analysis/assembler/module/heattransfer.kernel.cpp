
#include "heattransfer.h"
#include "analysis/assembler/module/assembler.hpp"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "analysis/scheme/steadystate.h"
#include "math/physics/matrix_distributed.h"

#include <numeric>
#include <algorithm>

namespace espreso {

template <size_t gps, size_t ndim, enum ThermalConductivityConfiguration::MODEL model, bool direct, class Physics> struct SetConductivity;

template <size_t gps, size_t ndim, bool direct, class Physics> struct SetConductivity<gps, ndim, ThermalConductivityConfiguration::MODEL::ISOTROPIC, direct, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		if (direct) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][0][s] = value; }));
		} else {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][0][s] = value; }));
		}
	}
};

template <size_t gps, bool direct, class Physics> struct SetConductivity<gps, 2, ThermalConductivityConfiguration::MODEL::DIAGONAL, direct, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		if (direct) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][1][s] = value; }));
		} else {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][1][s] = value; }));
		}
	}
};

template <size_t gps, bool direct, class Physics> struct SetConductivity<gps, 3, ThermalConductivityConfiguration::MODEL::DIAGONAL, direct, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		Evaluator *kzz = subkernels.conductivity.conductivity->values.get(2, 2).evaluator;
		if (direct) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][2][s] = value; }));
		} else {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][2][s] = value; }));
		}
	}
};

template <size_t gps, bool direct, class Physics> struct SetConductivity<gps, 2, ThermalConductivityConfiguration::MODEL::SYMMETRIC, direct, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kxy = subkernels.conductivity.conductivity->values.get(0, 1).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		if (direct) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][2][s] = value; }));
		} else {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][2][s] = value; }));
		}
	}
};

template <size_t gps, bool direct, class Physics> struct SetConductivity<gps, 3, ThermalConductivityConfiguration::MODEL::SYMMETRIC, direct, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kxy = subkernels.conductivity.conductivity->values.get(0, 1).evaluator;
		Evaluator *kxz = subkernels.conductivity.conductivity->values.get(0, 2).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		Evaluator *kyz = subkernels.conductivity.conductivity->values.get(1, 2).evaluator;
		Evaluator *kzz = subkernels.conductivity.conductivity->values.get(2, 2).evaluator;
		if (direct) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][2][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][3][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][4][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][5][s] = value; }));
		} else {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][2][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][3][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][4][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][5][s] = value; }));
		}
	}
};

template <size_t gps, bool direct, class Physics> struct SetConductivity<gps, 2, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, direct, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kxy = subkernels.conductivity.conductivity->values.get(0, 1).evaluator;
		Evaluator *kyx = subkernels.conductivity.conductivity->values.get(1, 0).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		if (direct) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][2][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][3][s] = value; }));
		} else {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][2][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][3][s] = value; }));
		}
	}
};

template <size_t gps, bool direct, class Physics> struct SetConductivity<gps, 3, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, direct, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		Evaluator *kxx = subkernels.conductivity.conductivity->values.get(0, 0).evaluator;
		Evaluator *kxy = subkernels.conductivity.conductivity->values.get(0, 1).evaluator;
		Evaluator *kxz = subkernels.conductivity.conductivity->values.get(0, 2).evaluator;
		Evaluator *kyx = subkernels.conductivity.conductivity->values.get(1, 0).evaluator;
		Evaluator *kyy = subkernels.conductivity.conductivity->values.get(1, 1).evaluator;
		Evaluator *kyz = subkernels.conductivity.conductivity->values.get(1, 2).evaluator;
		Evaluator *kzx = subkernels.conductivity.conductivity->values.get(2, 0).evaluator;
		Evaluator *kzy = subkernels.conductivity.conductivity->values.get(2, 1).evaluator;
		Evaluator *kzz = subkernels.conductivity.conductivity->values.get(2, 2).evaluator;
		if (direct) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][2][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][3][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][4][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][5][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][6][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][7][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.conductivity[gp][8][s] = value; }));
		} else {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kxz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][2][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][3][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][4][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kyz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][5][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzx, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][6][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzy, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][7][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				kzz, [] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.conductivity[gp][8][s] = value; }));
		}
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetTranslation;

template <size_t gps, class Physics> struct SetTranslation<gps, 2, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		if (subkernels.coosystem.configuration) {
			switch (subkernels.coosystem.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				break;
			}
		}
	}
};

template <size_t gps, class Physics> struct SetTranslation<gps, 3, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		if (subkernels.coosystem.configuration) {
			switch (subkernels.coosystem.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->rotation.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.x.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][0][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.y.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][1][s] = value; }));
				subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
						subkernels.coosystem.configuration->center.z.evaluator,
						[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.center[gp][2][s] = value; }));
				break;
			}
		}
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetThickness {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{

	}
};

template <size_t gps, class Physics> struct SetThickness<gps, 2, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		if (subkernels.thickness.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.thickness.expression->evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.thickness[gp][s] = value; }));
		}
	}
};

template <size_t gps, size_t ndim, class Physics> struct SetAdvection;

template <size_t gps, class Physics> struct SetAdvection<gps, 2, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		if (subkernels.advection.isactive) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.advection.expression->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.advection.expression->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][1][s] = value; }));
		}
	}
};

template <size_t gps, class Physics> struct SetAdvection<gps, 3, Physics> {
	static void analyze(HeatTransfer::SubKernels &subkernels)
	{
		if (subkernels.advection.expression) {
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.advection.expression->x.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][0][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.advection.expression->y.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][1][s] = value; }));
			subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
					subkernels.advection.expression->z.evaluator,
					[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.advection[gp][2][s] = value; }));
		}
	}
};

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model>
void preprocess(HeatTransfer::SubKernels &subkernels)
{
	typedef HeatTransferElementDescriptor<nodes, gps, ndim, edim, ecfmodel, model> Physics;
	SetThickness<gps, ndim, Physics>::analyze(subkernels);
	if (subkernels.conductivity.indirect) {
		SetConductivity<gps, ndim, ecfmodel, false, Physics>::analyze(subkernels);
	} else {
		SetConductivity<gps, ndim, model, true, Physics>::analyze(subkernels);
	}
	if (subkernels.coosystem.isactive) {
		SetTranslation<gps, ndim, Physics>::analyze(subkernels);
	}
	if (subkernels.advection.isactive) {
		SetAdvection<gps, ndim, Physics>::analyze(subkernels);
	}

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);

	typename Physics::Element element;
	basis.simd(element);
	SIMD volume;
	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		integration.simd(element);
		for (size_t gp = 0; gp < gps; ++gp) {
			volume = volume + element.det[gp] * load1(element.w[gp]);
		}
	}

	subkernels.esize = sizeof(typename Physics::Element);
	for (size_t s = 0; s < SIMD::size; ++s) {
		subkernels.volume += volume[s];
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model>
void compute(const HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	typedef HeatTransferElementDescriptor<nodes, gps, ndim, edim, ecfmodel, model> Physics;
	typename Physics::Element element;

	BasisKernel<code, nodes, gps, edim> basis(subkernels.basis);
	CoordinatesKernel<nodes, gps, ndim, Physics> coordinates(subkernels.coordinates);
	ThicknessToNodes<nodes, ndim, Physics> thickness(subkernels.thickness);
	TemperatureKernel<nodes, gps, Physics> temperature(subkernels.temperature);
	IntegrationKernel<nodes, gps, ndim, edim, Physics> integration(subkernels.integration);
	HeatTransferCoordinateSystemKernel<gps, ndim, ecfmodel, model, Physics> coosystem(subkernels.coosystem);
//	AdvectionKernel<nodes, gps, ndim, etype, Physics> advection(subkernels[interval].advection);
	HeatTransferMatrixKernel<nodes, gps, ndim, model, Physics> K(subkernels.K);
	TemperatureGradientKernel<nodes, gps, ndim, Physics> gradient(subkernels.gradient);
	TemperatureFluxKernel<nodes, gps, ndim, model, Physics> flux(subkernels.flux);

	std::vector<ExternalGPsExpression<gps, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalGPsExpression<gps, Physics>*>(subkernels.expressions[i]));
		}
	}

	basis.simd(element);
	if (coosystem.isactive) {
		coosystem.simd(element);
	}

	thickness.setActiveness(action);
	temperature.setActiveness(action);
	coosystem.setActiveness(action);
//	advection.setActiveness(action);
	K.setActiveness(action);
	gradient.setActiveness(action);
	flux.setActiveness(action);

//	printf("sub-kernels: ");
	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
//		if (c == 0) printf("coordinates ");
		if (temperature.isactive) {
			temperature.simd(element);
//			if (c == 0) printf("temp ");
		}
		if (thickness.isactive) {
			thickness.simd(element);
//			if (c == 0) printf("thickness ");
		}
		integration.simd(element);
//		if (c == 0) printf("integrate ");
		if (coosystem.isactive) {
			coosystem.simd(element);
//			if (c == 0) printf("coosystem ");
		}
//		if (advection.isactive) {
//			advection.simd(element);
//		}
		if (K.isactive) {
			K.simd(element);
//			if (c == 0) printf("K ");
		}
		if (gradient.isactive) {
			gradient.simd(element);
//			if (c == 0) printf("gradient ");
		}
		if (flux.isactive) {
			flux.simd(element);
//			if (c == 0) printf("flux ");
		}
	}
//	printf("\n");
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model>
void fill(const HeatTransfer::SubKernels &subkernels)
{
	typedef HeatTransferElementDescriptor<nodes, gps, ndim, edim, ecfmodel, model> Physics;
	typename Physics::Element element;

	MatricFillerKernel<nodes, Physics> K(subkernels.Kfiller);
	VectorFillerKernel<nodes, Physics> RHS(subkernels.RHSfiller);

	for (esint c = 0; c < subkernels.chunks; ++c) {
		if (K.isactive) {
			K.simd(element);
		}
		if (RHS.isactive) {
			RHS.simd(element);
		}
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void preprocess(HeatTransfer::BoundarySubKernels &subkernels)
{

}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void compute(const HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void fill(const HeatTransfer::BoundarySubKernels &subkernels)
{

}

template <size_t ndim>
void initDirichlet(HeatTransfer::BoundarySubKernels &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<1, 1, ndim, 0> Physics;
	if (subkernels.temperature.expression) {
		subkernels.expressions.push_back(new ExternalNodeExpression<1, Physics>(
				subkernels.temperature.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.temp[0][s] = value; }));
	}
}

template <size_t ndim>
void dirichlet(const HeatTransfer::BoundarySubKernels &subkernels)
{
	typedef HeatTransferBoundaryDescriptor<1, 1, ndim, 0> Physics;
	typename Physics::Element element;

	CoordinatesKernel<1, 1, ndim, Physics> coordinates(subkernels.coordinates);
	VectorSetterKernel<1, Physics> set(subkernels.dirichlet, [] (auto &element, size_t &n, size_t &d, size_t &s) { return element.temp[n][s]; });

	std::vector<ExternalNodeExpression<1, Physics>*> nonconst;
	for (size_t i = 0; i < subkernels.expressions.size(); ++i) {
		if (subkernels.expressions[i]->evaluator->isConst()) {
			dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i])->simd(element);
		} else {
			nonconst.push_back(dynamic_cast<ExternalNodeExpression<1, Physics>*>(subkernels.expressions[i]));
		}
	}

	for (esint c = 0; c < subkernels.chunks; ++c) {
		coordinates.simd(element);
		for (size_t i = 0; i < nonconst.size(); ++i) {
			nonconst[i]->simd(element);
		}
		set.simd(element);
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model>
void runAction(HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim, ecfmodel, model>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim, ecfmodel, model>(subkernels); break;
	default: compute<code, nodes, gps, ndim, edim, ecfmodel, model>(subkernels, action); break;
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runAction(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim>(subkernels); break;
	default: compute<code, nodes, gps, ndim, edim>(subkernels, action); break;
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runConductivity(HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.conductivity.conductivity->model) {
	case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		if (subkernels.advection.isactive) {
			runAction<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ISOTROPIC, ThermalConductivityConfiguration::MODEL::DIAGONAL>(subkernels, action);
		} else {
			runAction<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ISOTROPIC, ThermalConductivityConfiguration::MODEL::ISOTROPIC>(subkernels, action);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::DIAGONAL:
		if (subkernels.coosystem.rotated) {
			runAction<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::SYMMETRIC>(subkernels, action);
		} else {
			runAction<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::DIAGONAL, ThermalConductivityConfiguration::MODEL::DIAGONAL>(subkernels, action);
		}
		break;
	case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
		runAction<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::SYMMETRIC, ThermalConductivityConfiguration::MODEL::SYMMETRIC>(subkernels, action);
		break;
	case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
		runAction<code, nodes, gps, ndim, edim, ThermalConductivityConfiguration::MODEL::ANISOTROPIC, ThermalConductivityConfiguration::MODEL::ANISOTROPIC>(subkernels, action);
		break;
	}
}

template <size_t ndim> void runElement(HeatTransfer::SubKernels &subkernels, Assembler::Action action);

template <> void runElement<2>(HeatTransfer::SubKernels &subkernels, Assembler::Action action)

{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): runConductivity<Element::CODE::TRIANGLE3, 3, HeatTransferGPC::TRIANGLE3, 2, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): runConductivity<Element::CODE::TRIANGLE6, 6, HeatTransferGPC::TRIANGLE6, 2, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE4  ): runConductivity<Element::CODE::SQUARE4  , 4, HeatTransferGPC::SQUARE4  , 2, 2>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::SQUARE8  ): runConductivity<Element::CODE::SQUARE8  , 8, HeatTransferGPC::SQUARE8  , 2, 2>(subkernels, action); break;
	}
}

template <> void runElement<3>(HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	switch (subkernels.code) {
	case static_cast<size_t>(Element::CODE::TETRA4   ): runConductivity<Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4    , 3, 3>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::TETRA10  ): runConductivity<Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10   , 3, 3>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5 ): runConductivity<Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5  , 3, 3>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): runConductivity<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13 , 3, 3>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PRISMA6  ): runConductivity<Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6   , 3, 3>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::PRISMA15 ): runConductivity<Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15  , 3, 3>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::HEXA8    ): runConductivity<Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8     , 3, 3>(subkernels, action); break;
	case static_cast<size_t>(Element::CODE::HEXA20   ): runConductivity<Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20    , 3, 3>(subkernels, action); break;
	}
}

template <size_t ndim, size_t edim> void runBoundary(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action);

template <> void runBoundary<2, 1>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <> void runBoundary<2, 0>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<2>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<2>(subkernels); break;
	}
}
template <> void runBoundary<3, 2>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <> void runBoundary<3, 1>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{

}

template <> void runBoundary<3, 0>(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: initDirichlet<3>(subkernels); break;
	case Assembler::Action::ASSEMBLE: case Assembler::Action::REASSEMBLE: dirichlet<3>(subkernels); break;
	}
}

void HeatTransfer::run(Action action, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2: runElement<2>(subkernels[interval], action); break;
	case 3: runElement<3>(subkernels[interval], action); break;
	}
}

void HeatTransfer::run(Action action, size_t region, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: runBoundary<2, 0>(boundary[region][interval], action); break;
		case 1: runBoundary<2, 1>(boundary[region][interval], action); break;
		} break;
	case 3:
		switch (info::mesh->boundaryRegions[region]->dimension) {
		case 0: runBoundary<3, 0>(boundary[region][interval], action); break;
		case 1: runBoundary<3, 1>(boundary[region][interval], action); break;
		case 2: runBoundary<3, 2>(boundary[region][interval], action); break;
		} break;
	}
}

}
