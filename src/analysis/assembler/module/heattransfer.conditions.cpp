
#include "heattransfer.h"
#include "heattransfer.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/temperature.h"
#include "analysis/assembler/operators/advection.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/conductivity.coordinatesystem.h"
#include "analysis/assembler/operators/heattransfer.f.h"
#include "analysis/assembler/operators/heattransfer.K.h"
#include "analysis/assembler/operators/filler.h"
#include "analysis/assembler/operators/gradient.h"
#include "analysis/assembler/operators/flux.h"

#include "basis/expression/variable.h"
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateConductivity {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateConductivity<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.conductivity[gp][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateConductivity<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->thermal_conductivity.model) {
		case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
			break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:
			mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			mat->thermal_conductivity.values.get(0, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
				}
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateConductivity<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->thermal_conductivity.model) {
		case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
			break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:
			mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][4][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(2, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][8][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			mat->thermal_conductivity.values.get(0, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(0, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][6][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][5][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][7][s] = results[gps * s + gp];
				}
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateConductivity<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.conductivity[gp][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateConductivity<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->thermal_conductivity.model) {
		case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
			mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(0, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:
			mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			mat->thermal_conductivity.values.get(0, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
				}
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateConductivity<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->thermal_conductivity.model) {
		case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
		case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
			mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][4][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(2, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][8][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(0, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(0, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][5][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(2, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][6][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(2, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][7][s] = results[gps * s + gp];
				}
			}
			break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:
			mat->thermal_conductivity.values.get(0, 0).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][4][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(2, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][8][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			mat->thermal_conductivity.values.get(0, 1).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(0, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][6][s] = results[gps * s + gp];
				}
			}
			mat->thermal_conductivity.values.get(1, 2).evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][5][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][7][s] = results[gps * s + gp];
				}
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateRotation {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateRotation<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			mat->coordinate_system.rotation.z.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			mat->coordinate_system.center.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateRotation<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			mat->coordinate_system.rotation.z.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			mat->coordinate_system.center.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateRotation<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			mat->coordinate_system.rotation.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.rotation.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.rotation.z.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			mat->coordinate_system.center.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			mat->coordinate_system.center.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.z.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateRotation<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN:
			mat->coordinate_system.rotation.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.rotation.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.rotation.z.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			mat->coordinate_system.center.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			mat->coordinate_system.center.x.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.y.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			mat->coordinate_system.center.z.evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct applyRotation {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {

	applyRotation(size_t interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {

	applyRotation(size_t interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCylindric;
	HeatTransferCoordinateSystemSpherical<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationSpherical;
	HeatTransferCoordinateSystemApply    <nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > apply;

	applyRotation(size_t interval): rotationCartesian(interval), rotationCylindric(interval), rotationSpherical(interval), apply(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL  : rotationSpherical.simd(element); break;
		}
		apply.simd(element);
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCylindric;
	HeatTransferCoordinateSystemApply    <nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > apply;

	applyRotation(size_t interval): rotationCartesian(interval), rotationCylindric(interval), apply(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		}
		apply.simd(element);
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCylindric;
	HeatTransferCoordinateSystemSpherical<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationSpherical;

	applyRotation(size_t interval): rotationCartesian(interval), rotationCylindric(interval), rotationSpherical(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL  : rotationSpherical.simd(element); break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCylindric;

	applyRotation(size_t interval): rotationCartesian(interval), rotationCylindric(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateHeatSource {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.heatSource[gp][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateTranslationMotion {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, Evaluator* evaluator)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateTranslationMotion<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.advection[gp][0][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateTranslationMotion<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		evaluator->evalVector(SIMD::size * gps, Evaluator::Params(), results);
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.advection[gp][0][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
Assembler::measurements HeatTransfer::conditionsloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	double initStart, initEnd;
	initStart = eslog::time();
	if (this->K == nullptr) {
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, etype>(action, ops, elements);
	}
	if (elements == 0) return {0.0, 0.0};
	typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element element;

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if (elements > SIMD::size) {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->simd(element);
			} else {
				dynamic_cast<DataDescriptor<nodes, gps, ndim, edim, etype>*>(*op)->peel(element, elements);
			}
		}
	}

	auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[interval].begin;
	CoordinatesToElementNodes<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > coo(interval, procNodes);
	CoordinatesToElementNodesAndGPs<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > cooAndGps(interval, procNodes);
	TemperatureToElementNodes<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > temp(interval, procNodes, Results::temperature->data.data());
	Integration<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > integration(interval);
	HeatTransferStiffness<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > stiffness(interval, this->elements.stiffness);
	HeatSource<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > heatSource(interval, this->elements.rhs);
	Advection<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > advection(interval, this->elements.stiffness);

	applyRotation<DataDescriptor, nodes, gps, ndim, edim, etype> rotation(interval);
	updateHeatSource<DataDescriptor, nodes, gps, ndim, edim, etype> updateHS;
	updateTranslationMotion<DataDescriptor, nodes, gps, ndim, edim, etype> updateTM;

	SymmetricMatricFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > upperFiller(interval, 1, this->elements.stiffness, this->K);
	GeneralMatricFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > fullFiller(interval, 1, this->elements.stiffness, this->K);
	VectorFiller<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > rhsFiller(interval, 1, this->elements.rhs, this->f);

	TemperatureGradient<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > gradient(interval, Results::gradient);
	TemperatureFlux<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > flux(interval, Results::flux);

	coo.move(SIMD::size);
	cooAndGps.move(SIMD::size);
	temp.move(SIMD::size);
	stiffness.move(SIMD::size);
	advection.move(SIMD::size);
	upperFiller.move(SIMD::size);
	fullFiller.move(SIMD::size);
	rhsFiller.move(SIMD::size);
	heatSource.move(SIMD::size);

	gradient.move(ndim * SIMD::size);
	flux.move(ndim * SIMD::size);

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator->params.general.size() == 0;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator->params.general.size() == 0 && advectionEval->second.y.evaluator->params.general.size() == 0 && advectionEval->second.z.evaluator->params.general.size() == 0;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;
	bool isfullMatrix = this->K->shape == Matrix_Shape::FULL;

	initEnd = eslog::time();
	double start, end;

	if(action == ActionOperator::ASSEMBLE)
	{
		start = eslog::time();
		__SSC_MARK(0xFACE);
		esint chunks = elements / SIMD::size;
		for (esint c = 1; c < chunks; ++c) {
			if (cooToGP) {
				cooAndGps.simd(element);
			} else {
				coo.simd(element);
			}
			integration.simd(element);
			if (getTemp) {
				temp.simd(element);
			}
			if (computeConductivity) {
				if (!constConductivity) {
					updateConductivity<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
				}
				if (rotateConductivity) {
					if (!constRotation) {
						updateRotation<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
					}
					rotation(element, mat);
				}
			}

			if (computeK) {
				if (hasAdvection) {
					if (!constAdvection) {
						updateTM(element, advectionEval->second.x.evaluator);
					}
					advection.simd(element);
				}
				stiffness.simd(element);
				if (hasHeatSource) {
					if (!constHeatSource) {
						updateHS(element, heatSourceEval->second.evaluator);
					}
					heatSource.simd(element);
				}
			}
			if (action == ActionOperator::FILL) {
				if (isfullMatrix) {
					fullFiller.simd(element);
				} else {
					upperFiller.simd(element);
				}
				rhsFiller.simd(element);
			}
			if (computeGradient) {
				gradient.simd(element);
			}
			if (computeFlux) {
				flux.simd(element);
			}
		}
		__SSC_MARK(0xDEAD);
		end = eslog::time();
	}

	if(action == ActionOperator::REASSEMBLE)
	{
		start = eslog::time();
		__SSC_MARK(0xCAFE);
		esint chunks = elements / SIMD::size;
		for (esint c = 1; c < chunks; ++c) {
			if (cooToGP) {
				cooAndGps.simd(element);
			} else {
				coo.simd(element);
			}
			integration.simd(element);
			if (getTemp) {
				temp.simd(element);
			}
			if (computeConductivity) {
				if (!constConductivity) {
					updateConductivity<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
				}
				if (rotateConductivity) {
					if (!constRotation) {
						updateRotation<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
					}
					rotation(element, mat);
				}
			}

			if (computeK) {
				if (hasAdvection) {
					if (!constAdvection) {
						updateTM(element, advectionEval->second.x.evaluator);
					}
					advection.simd(element);
				}
				stiffness.simd(element);
				if (hasHeatSource) {
					if (!constHeatSource) {
						updateHS(element, heatSourceEval->second.evaluator);
					}
					heatSource.simd(element);
				}
			}
			if (action == ActionOperator::FILL) {
				if (isfullMatrix) {
					fullFiller.simd(element);
				} else {
					upperFiller.simd(element);
				}
				rhsFiller.simd(element);
			}
			if (computeGradient) {
				gradient.simd(element);
			}
			if (computeFlux) {
				flux.simd(element);
			}
		}
		__SSC_MARK(0xFADE);
		end = eslog::time();
	}

	if((action != ActionOperator::REASSEMBLE) && (action != ActionOperator::ASSEMBLE))
	{
		start = eslog::time();
		esint chunks = elements / SIMD::size;
		for (esint c = 1; c < chunks; ++c) {
			if (cooToGP) {
				cooAndGps.simd(element);
			} else {
				coo.simd(element);
			}
			integration.simd(element);
			if (getTemp) {
				temp.simd(element);
			}
			if (computeConductivity) {
				if (!constConductivity) {
					updateConductivity<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
				}
				if (rotateConductivity) {
					if (!constRotation) {
						updateRotation<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
					}
					rotation(element, mat);
				}
			}

			if (computeK) {
				if (hasAdvection) {
					if (!constAdvection) {
						updateTM(element, advectionEval->second.x.evaluator);
					}
					advection.simd(element);
				}
				stiffness.simd(element);
				if (hasHeatSource) {
					if (!constHeatSource) {
						updateHS(element, heatSourceEval->second.evaluator);
					}
					heatSource.simd(element);
				}
			}
			if (action == ActionOperator::FILL) {
				if (isfullMatrix) {
					fullFiller.simd(element);
				} else {
					upperFiller.simd(element);
				}
				rhsFiller.simd(element);
			}
			if (computeGradient) {
				gradient.simd(element);
			}
			if (computeFlux) {
				flux.simd(element);
			}
		}
		end = eslog::time();
	}

	if (elements % SIMD::size) {
		eslog::error("peel loop is not supported\n");
		// peel is never needed
	}

	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if ((*op)->isconst) {
				(*op)->move(-(int)std::min(elements, (esint)SIMD::size));
			} else {
				(*op)->move(-(esint)SIMD::size);
			}
		}
	}
	return {initEnd - initStart,  end - start};
}

template <int etype>
Assembler::measurements HeatTransfer::instantiateConditions2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return conditionsloop<HeatTransferDataDescriptor, 3, HeatTransferGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return conditionsloop<HeatTransferDataDescriptor, 6, HeatTransferGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return conditionsloop<HeatTransferDataDescriptor, 4, HeatTransferGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return conditionsloop<HeatTransferDataDescriptor, 8, HeatTransferGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
	default: return { .0, .0 };
	};
}

template <int etype>
Assembler::measurements HeatTransfer::instantiateConditions3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return conditionsloop<HeatTransferDataDescriptor,  4, HeatTransferGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return conditionsloop<HeatTransferDataDescriptor, 10, HeatTransferGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return conditionsloop<HeatTransferDataDescriptor,  5, HeatTransferGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return conditionsloop<HeatTransferDataDescriptor, 13, HeatTransferGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return conditionsloop<HeatTransferDataDescriptor,  6, HeatTransferGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return conditionsloop<HeatTransferDataDescriptor, 15, HeatTransferGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return conditionsloop<HeatTransferDataDescriptor,  8, HeatTransferGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return conditionsloop<HeatTransferDataDescriptor, 20, HeatTransferGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
	default: return { .0, .0 };
	};
}

Assembler::measurements HeatTransfer::instantiateConditions(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiateConditions2D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiateConditions2D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiateConditions2D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiateConditions2D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case HeatTransferElementType::SYMMETRIC_ISOTROPIC : return instantiateConditions3D<HeatTransferElementType::SYMMETRIC_ISOTROPIC >(action, code, ops, interval, elements);
		case HeatTransferElementType::SYMMETRIC_GENERAL   : return instantiateConditions3D<HeatTransferElementType::SYMMETRIC_GENERAL   >(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_ISOTROPIC: return instantiateConditions3D<HeatTransferElementType::ASYMMETRIC_ISOTROPIC>(action, code, ops, interval, elements);
		case HeatTransferElementType::ASYMMETRIC_GENERAL  : return instantiateConditions3D<HeatTransferElementType::ASYMMETRIC_GENERAL  >(action, code, ops, interval, elements);
		}
	}
	return {0.0, 0.0};
}

}
