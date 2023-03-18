
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
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
		}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate();
			}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][4][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][8][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][6][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 2).evaluator->evaluate();
			}
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
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
		}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate();
			}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][4][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][8][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][5][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(2, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][6][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(2, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][7][s] = results[gps * s + gp];
				}
			}
			break;
		case ThermalConductivityConfiguration::MODEL::DIAGONAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][4][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(2, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][8][s] = results[gps * s + gp];
				}
			}
			/* no break */
		case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][1][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][3][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(0, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.conductivity[gp][2][s] = results[gps * s + gp];
					element.ecf.conductivity[gp][6][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->thermal_conductivity.values.get(1, 2).evaluator->evaluate();
			}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.z.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.y.evaluator->evaluate();
			}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.z.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.y.evaluator->evaluate();
			}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.y.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.z.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.y.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.y.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.z.evaluator->evaluate();
			}
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
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.y.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.rotation.z.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.y.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.x.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.y.evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.center[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->coordinate_system.center.z.evaluator->evaluate();
			}
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
struct updateCosSin {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateCosSin<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {

	updateCosSin(size_t interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct updateCosSin<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {

	updateCosSin(size_t interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateCosSin<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCylindric;
	HeatTransferCoordinateSystemSpherical<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationSpherical;

	updateCosSin(size_t interval): rotationCartesian(interval), rotationCylindric(interval), rotationSpherical(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL  : rotationSpherical.simd(element); break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateCosSin<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > rotationCylindric;

	updateCosSin(size_t interval): rotationCartesian(interval), rotationCylindric(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateCosSin<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCylindric;
	HeatTransferCoordinateSystemSpherical<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationSpherical;

	updateCosSin(size_t interval): rotationCartesian(interval), rotationCylindric(interval), rotationSpherical(interval) {}

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
struct updateCosSin<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemCartesian<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCartesian;
	HeatTransferCoordinateSystemCylindric<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > rotationCylindric;

	updateCosSin(size_t interval): rotationCartesian(interval), rotationCylindric(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct applyRotation {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC> {

	applyRotation(size_t interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_ISOTROPIC>::Element &element)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC> {

	applyRotation(size_t interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_ISOTROPIC>::Element &element)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemApply <nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_GENERAL, DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_GENERAL> > apply;

	applyRotation(size_t interval): apply(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element)
	{
		apply.simd(element);
	}
};


template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim>
struct applyRotation<DataDescriptor, nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {
	HeatTransferCoordinateSystemApply <nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_GENERAL, DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> > apply;

	applyRotation(size_t interval): apply(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element)
	{
		apply.simd(element);
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateHeatSource {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = evaluator->evaluate();
		}
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
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = evaluator->evaluate();
		}
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
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = evaluator->evaluate();
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.advection[gp][0][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct storeCosSin {

	storeCosSin(std::vector<double> &storage, size_t elements) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct storeCosSin<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	storeCosSin(std::vector<double> &storage, size_t elements): storage(storage)
	{
		storage.resize(2 * elements * gps * SIMD::size);
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(iterator, element.cossin, 2 * gps * SIMD::size * sizeof(double));
		iterator += 2 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct storeCosSin<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	storeCosSin(std::vector<double> &storage, size_t elements): storage(storage)
	{
		storage.resize(6 * elements * gps * SIMD::size);
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(iterator, element.cossin, 6 * gps * SIMD::size * sizeof(double));
		iterator += 6 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct storeCosSin<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	storeCosSin(std::vector<double> &storage, size_t elements): storage(storage)
	{
		storage.resize(2 * elements * gps * SIMD::size);
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(iterator, element.cossin, 2 * gps * SIMD::size * sizeof(double));
		iterator += 2 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct storeCosSin<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	storeCosSin(std::vector<double> &storage, size_t elements): storage(storage)
	{
		storage.resize(6 * elements * gps * SIMD::size);
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(iterator, element.cossin, 6 * gps * SIMD::size * sizeof(double));
		iterator += 6 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct loadCosSin {

	loadCosSin(std::vector<double> &storage) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct loadCosSin<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	loadCosSin(std::vector<double> &storage): storage(storage)
	{
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(element.cossin, iterator, 2 * gps * SIMD::size * sizeof(double));
		iterator += 2 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct loadCosSin<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	loadCosSin(std::vector<double> &storage): storage(storage)
	{
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::SYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(element.cossin, iterator, 6 * gps * SIMD::size * sizeof(double));
		iterator += 6 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct loadCosSin<DataDescriptor, nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	loadCosSin(std::vector<double> &storage): storage(storage)
	{
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(element.cossin, iterator, 2 * gps * SIMD::size * sizeof(double));
		iterator += 2 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct loadCosSin<DataDescriptor, nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL> {

	std::vector<double> &storage;
	double* iterator;
	loadCosSin(std::vector<double> &storage): storage(storage)
	{
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, HeatTransferElementType::ASYMMETRIC_GENERAL>::Element &element)
	{
		memcpy(element.cossin, iterator, 6 * gps * SIMD::size * sizeof(double));
		iterator += 6 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateThickness {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, Evaluator* evaluator)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateThickness<DataDescriptor, nodes, gps, 2, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = evaluator->evaluate();
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.thickness[gp][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateThickness<DataDescriptor, nodes, gps, 3, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element, Evaluator* evaluator)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
Assembler::measurements HeatTransfer::conditionsloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	if (this->K == nullptr) {
		return loop<HeatTransferDataDescriptor, nodes, gps, ndim, edim, etype>(action, ops, elements);
	}
	if (elements == 0) return {0.0, 0.0};

	double initStart = eslog::time();

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

	updateCosSin<DataDescriptor, nodes, gps, ndim, edim, etype> cossin(interval);
	storeCosSin<DataDescriptor, nodes, gps, ndim, edim, etype> storecossin(cossin_conditions[interval], elements);
	loadCosSin<DataDescriptor, nodes, gps, ndim, edim, etype> loadcossin(cossin_conditions[interval]);
	applyRotation<DataDescriptor, nodes, gps, ndim, edim, etype> rotation(interval);
	updateHeatSource<DataDescriptor, nodes, gps, ndim, edim, etype> updateHS;
	updateTranslationMotion<DataDescriptor, nodes, gps, ndim, edim, etype> updateTM;
	updateThickness<DataDescriptor, nodes, gps, ndim, edim, etype> updateThickness;

	TemperatureGradient<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > gradient(interval, Results::gradient);
	TemperatureFlux<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > flux(interval, Results::flux);

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateConductivity = mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
	bool constCosSin = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constConductivity = !settings.diffusion_split;
	bool constRotation = settings.sigma == 0;
	if (mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateConductivity &= mat->coordinate_system.rotation.z.evaluator->parameters.size();
			}
			if (ndim == 3) {
				rotateConductivity &= mat->coordinate_system.rotation.x.evaluator->parameters.size() | mat->coordinate_system.rotation.y.evaluator->parameters.size() | mat->coordinate_system.rotation.z.evaluator->parameters.size();
			}
		}
	}
	if (info::ecf->always_update_conductivity) { // TODO
		constConductivity = false;
	}
	bool storeCosSin = settings.reassembling_optimization && action == ActionOperator::ASSEMBLE;
	bool loadCosSin  = settings.reassembling_optimization && action != ActionOperator::ASSEMBLE;

	bool hasHeatSource = false;
	bool constHeatSource = true;
	auto heatSourceEval = configuration.heat_source.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (heatSourceEval != configuration.heat_source.end()) {
		hasHeatSource = true;
		constHeatSource = heatSourceEval->second.evaluator != nullptr;
	}

	bool hasAdvection = false;
	bool constAdvection = true;
	auto advectionEval = configuration.translation_motions.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (advectionEval != configuration.translation_motions.end()) {
		hasAdvection = true;
		constAdvection = advectionEval->second.x.evaluator != nullptr && advectionEval->second.y.evaluator != nullptr && advectionEval->second.z.evaluator != nullptr;
	}

	bool constThickness = true;
	auto thicknessEval = settings.thickness.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (thicknessEval != settings.thickness.end()) {
		constThickness = thicknessEval->second.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeGradient = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.gradient;
	bool computeFlux = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.flux;
	bool computeConductivity = computeK | computeFlux;
	bool getTemp = computeGradient || computeFlux;

	esint chunks = elements / SIMD::size;
	double init = eslog::time() - initStart;

	if (cooToGP) {
		cooAndGps.move(SIMD::size);
	} else {
		coo.move(SIMD::size);
	}
	if (getTemp) {
		temp.move(SIMD::size);
	}
	if (computeK) {
		stiffness.move(SIMD::size);
		if (hasHeatSource) {
			heatSource.move(SIMD::size);
		}
	}
	if (hasAdvection) {
		advection.move(SIMD::size);
	}
	if (computeGradient) {
		gradient.move(SIMD::size);
	}
	if (computeFlux) {
		flux.move(SIMD::size);
	}

	double start = eslog::time();
	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xFACE); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xCAFE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xCAFE); break; // TODO
	default:
		eslog::error("unsupported action\n");
	}
	for (esint c = 1; c < chunks; ++c) {
		if (cooToGP) {
			cooAndGps.simd(element);
//			printf(" cooAndGps");
		} else {
			coo.simd(element);
//			printf(" coo");
		}
		if (!constThickness) {
			updateThickness(element, thicknessEval->second.evaluator);
//			printf(" thickness");
		}

		integration.simd(element);
//		printf(" integration");
		if (getTemp) {
			temp.simd(element);
//			printf(" temp");
		}
		if (computeConductivity) {
			if (!constConductivity) {
				updateConductivity<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
//				printf(" updateConductivity");
			}
			if (rotateConductivity) {
				if (!constRotation) {
					updateRotation<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
//					printf(" updateRotation");
				}
				if (!constCosSin) {
					if (loadCosSin) {
//						printf(" loadCosSin");
						loadcossin(element);
					} else {
//						printf(" updateCosSin");
						cossin(element, mat);
					}
					if (storeCosSin) {
//						printf(" storeCosSin");
						storecossin(element);
					}
				}
				rotation(element);
//				printf(" rotation");
			}
		}

		if (computeK) {
			if (hasAdvection) {
				if (!constAdvection) {
					updateTM(element, advectionEval->second.x.evaluator);
//					printf(" updateTM");
				}
				advection.simd(element);
//				printf(" advection");
			}
			stiffness.simd(element);
//			printf(" stiffness");
			if (hasHeatSource) {
				if (!constHeatSource) {
					updateHS(element, heatSourceEval->second.evaluator);
//					printf(" updateHS");
				}
				heatSource.simd(element);
//				printf(" heatSource");
			}
		}
		if (computeGradient) {
			gradient.simd(element);
//			printf(" gradient");
		}
		if (computeFlux) {
			flux.simd(element);
//			printf(" flux");
		}
//		printf("\n");
	}

	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xDEAD); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xFADE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xFADE); break; // TODO
	default:
		eslog::error("unsupported action\n");
	}
	double loop = eslog::time() - start;

	if (elements % SIMD::size) {
		eslog::error("peel loop is not supported\n");
		// peel is never needed
	}

	// move operators that was used for initialization only
	for (auto op = ops.cbegin(); op != ops.cend(); ++op) {
		if ((*op)->action & action) {
			if ((*op)->isconst) {
				(*op)->move(-std::min(elements, (esint)SIMD::size));
			} else {
				(*op)->move(-(esint)SIMD::size);
			}
		}
	}
	return { init, loop };
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
