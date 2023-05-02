
#include "analysis/assembler/module/heattransfer.h"
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

	if (subkernels.heatSource.isactive) {
		subkernels.expressions.push_back(new ExternalGPsExpression<gps, Physics>(
				subkernels.heatSource.expression->evaluator,
				[] (typename Physics::Element &element, size_t &gp, size_t &s, double value) { element.ecf.heatSource[gp][s] = value; }));
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
void fill(const HeatTransfer::BoundarySubKernels &subkernels)
{

}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim, enum ThermalConductivityConfiguration::MODEL ecfmodel, enum ThermalConductivityConfiguration::MODEL model>
void runAction(HeatTransfer::SubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim, ecfmodel, model>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim, ecfmodel, model>(subkernels); break;
	}
}

template <Element::CODE code, size_t nodes, size_t gps, size_t ndim, size_t edim>
void runAction(HeatTransfer::BoundarySubKernels &subkernels, Assembler::Action action)
{
	switch (action) {
	case Assembler::Action::PREPROCESS: preprocess<code, nodes, gps, ndim, edim>(subkernels); break;
	case Assembler::Action::FILL: fill<code, nodes, gps, ndim, edim>(subkernels); break;
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

void HeatTransfer::runPreprocess(Action action, size_t interval)
{
	switch (info::mesh->dimension) {
	case 2: runElement<2>(subkernels[interval], action); break;
	case 3: runElement<3>(subkernels[interval], action); break;
	}
}

}
