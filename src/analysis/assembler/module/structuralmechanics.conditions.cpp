
#include "structuralmechanics.h"
#include "structuralmechanics.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/displacement.h"
#include "analysis/assembler/operators/elasticity.h"
#include "analysis/assembler/operators/elasticity.coordinatesystem.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/structuralmechanics.f.h"
#include "analysis/assembler/operators/structuralmechanics.K.h"
#include "analysis/assembler/operators/stress.h"
#include "analysis/assembler/operators/filler.h"

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
struct updateElasticity {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateElasticity<DataDescriptor, nodes, gps, 2, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->linear_elastic_properties.model) {
		case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.youngModulus[gp][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.poissonRatio[gp][s] = results[gps * s + gp];
				}
			}
			break;
		case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
			break;
		case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
			break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim>
struct updateElasticity<DataDescriptor, nodes, gps, 3, edim, StructuralMechanicsElementType::SYMMETRIC_VOLUME> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, StructuralMechanicsElementType::SYMMETRIC_VOLUME>::Element &element, const MaterialConfiguration *mat)
	{
		double results[SIMD::size * gps];
		switch (mat->linear_elastic_properties.model) {
		case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.youngModulus[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.poissonRatio[gp][0][s] = results[gps * s + gp];
				}
			}
			break;
		case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.youngModulus[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.young_modulus.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.youngModulus[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.young_modulus.get(2, 2).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.youngModulus[gp][2][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.poissonRatio[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.poisson_ratio.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.poissonRatio[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.poisson_ratio.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.poissonRatio[gp][2][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.shear_modulus.get(0, 0).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.shearModulus[gp][0][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.shear_modulus.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.shearModulus[gp][1][s] = results[gps * s + gp];
				}
			}
			for (size_t v = 0; v < SIMD::size * gps; ++v) {
				results[v] = mat->linear_elastic_properties.shear_modulus.get(1, 1).evaluator->evaluate();
			}
			for (size_t gp = 0; gp < gps; ++gp) {
				for (size_t s = 0; s < SIMD::size; ++s) {
					element.ecf.shearModulus[gp][2][s] = results[gps * s + gp];
				}
			}
			break;
		case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
			break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateRotation {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateRotation<DataDescriptor, nodes, gps, 2, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element, const MaterialConfiguration *mat)
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateRotation<DataDescriptor, nodes, gps, 3, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element, const MaterialConfiguration *mat)
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateCosSin<DataDescriptor, nodes, gps, 3, edim, etype> {
	ElasticityCoordinateSystemCartesian<nodes, gps, 3, edim, etype, DataDescriptor<nodes, gps, 3, edim, etype> > rotationCartesian;
	ElasticityCoordinateSystemCylindric<nodes, gps, 3, edim, etype, DataDescriptor<nodes, gps, 3, edim, etype> > rotationCylindric;
	ElasticityCoordinateSystemSpherical<nodes, gps, 3, edim, etype, DataDescriptor<nodes, gps, 3, edim, etype> > rotationSpherical;

	updateCosSin(size_t interval): rotationCartesian(interval), rotationCylindric(interval), rotationSpherical(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element, const MaterialConfiguration *mat)
	{
		switch (mat->coordinate_system.type) {
		case CoordinateSystemConfiguration::TYPE::CARTESIAN  : rotationCartesian.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::CYLINDRICAL: rotationCylindric.simd(element); break;
		case CoordinateSystemConfiguration::TYPE::SPHERICAL  : rotationSpherical.simd(element); break;
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateCosSin<DataDescriptor, nodes, gps, 2, edim, etype> {
	ElasticityCoordinateSystemCartesian<nodes, gps, 2, edim, etype, DataDescriptor<nodes, gps, 2, edim, etype> > rotationCartesian;
	ElasticityCoordinateSystemCylindric<nodes, gps, 2, edim, etype, DataDescriptor<nodes, gps, 2, edim, etype> > rotationCylindric;

	updateCosSin(size_t interval): rotationCartesian(interval), rotationCylindric(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element, const MaterialConfiguration *mat)
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct applyRotation<DataDescriptor, nodes, gps, 2, edim, etype> {
	ElasticityCoordinateSystemApply <nodes, gps, 2, edim, etype, DataDescriptor<nodes, gps, 2, edim, etype> > apply;

	applyRotation(size_t interval): apply(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element)
	{
		apply.simd(element);
	}
};


template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct applyRotation<DataDescriptor, nodes, gps, 3, edim, etype> {
	ElasticityCoordinateSystemApply <nodes, gps, 3, edim, etype, DataDescriptor<nodes, gps, 3, edim, etype> > apply;

	applyRotation(size_t interval): apply(interval) {}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element)
	{
		apply.simd(element);
	}
};


template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct updateAcceleration {
	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = evaluator->evaluate();
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.acceleration[gp][0][s] = results[gps * s + gp];
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

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct storeCosSin<DataDescriptor, nodes, gps, 2, edim, etype> {

	std::vector<double> &storage;
	double* iterator;
	storeCosSin(std::vector<double> &storage, size_t elements): storage(storage)
	{
		storage.resize(4 * elements * gps * SIMD::size);
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element)
	{
		memcpy(iterator, element.cossin, 4 * gps * SIMD::size * sizeof(double));
		iterator += 4 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct storeCosSin<DataDescriptor, nodes, gps, 3, edim, etype> {

	std::vector<double> &storage;
	double* iterator;
	storeCosSin(std::vector<double> &storage, size_t elements): storage(storage)
	{
		storage.resize(12 * elements * gps * SIMD::size);
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element)
	{
		memcpy(iterator, element.cossin, 12 * gps * SIMD::size * sizeof(double));
		iterator += 12 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype>
struct loadCosSin {

	loadCosSin(std::vector<double> &storage) {}

	void operator()(typename DataDescriptor<nodes, gps, ndim, edim, etype>::Element &element)
	{

	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct loadCosSin<DataDescriptor, nodes, gps, 2, edim, etype> {

	std::vector<double> &storage;
	double* iterator;
	loadCosSin(std::vector<double> &storage): storage(storage)
	{
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element)
	{
		memcpy(element.cossin, iterator, 4 * gps * SIMD::size * sizeof(double));
		iterator += 4 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct loadCosSin<DataDescriptor, nodes, gps, 3, edim, etype> {

	std::vector<double> &storage;
	double* iterator;
	loadCosSin(std::vector<double> &storage): storage(storage)
	{
		iterator = storage.data();
	}

	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element)
	{
		memcpy(element.cossin, iterator, 12 * gps * SIMD::size * sizeof(double));
		iterator += 12 * gps * SIMD::size;
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t ndim, size_t edim, size_t etype> struct updateVelocity;

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateVelocity<DataDescriptor, nodes, gps, 2, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 2, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = evaluator->evaluate();
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.angularVelocity[gp][s] = results[gps * s + gp];
			}
		}
	}
};

template <template <size_t, size_t, size_t, size_t, size_t> class DataDescriptor, size_t nodes, size_t gps, size_t edim, size_t etype>
struct updateVelocity<DataDescriptor, nodes, gps, 3, edim, etype> {
	void operator()(typename DataDescriptor<nodes, gps, 3, edim, etype>::Element &element, Evaluator* evaluator)
	{
		double results[SIMD::size * gps];
		for (size_t v = 0; v < SIMD::size * gps; ++v) {
			results[v] = evaluator->evaluate();
		}
		for (size_t gp = 0; gp < gps; ++gp) {
			for (size_t s = 0; s < SIMD::size; ++s) {
				element.ecf.angularVelocity[gp][0][s] = results[gps * s + gp];
			}
		}
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
Assembler::measurements StructuralMechanics::conditionsloop(ActionOperator::Action action, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	if (this->K == nullptr) {
		return loop<StructuralMechanicsDataDescriptor, nodes, gps, ndim, edim, etype>(action, ops, elements);
	}
	if (elements == 0) return measurements();

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
	DisplacementToElementNodes<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > displacement(interval, procNodes, Results::displacement->data.data());
	Integration<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > integration(interval);
	StructuralMechanicsStiffness<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > stiffness(interval, this->elements.stiffness);
	Acceleration<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > acceleration(interval, this->elements.rhs);
	AngularVelocity<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > velocity(interval, this->elements.rhs);

	updateCosSin<DataDescriptor, nodes, gps, ndim, edim, etype> cossin(interval);
	storeCosSin<DataDescriptor, nodes, gps, ndim, edim, etype> storecossin(cossin_conditions[interval], elements);
	loadCosSin<DataDescriptor, nodes, gps, ndim, edim, etype> loadcossin(cossin_conditions[interval]);
	applyRotation<DataDescriptor, nodes, gps, ndim, edim, etype> rotation(interval);
	updateAcceleration<DataDescriptor, nodes, gps, ndim, edim, etype> updateAcc;
	updateVelocity<DataDescriptor, nodes, gps, ndim, edim, etype> updateVelocity;
	updateThickness<DataDescriptor, nodes, gps, ndim, edim, etype> updateThickness;

	Stress<nodes, gps, ndim, edim, etype, DataDescriptor<nodes, gps, ndim, edim, etype> > stress(interval, Results::principalStress, Results::componentStress, Results::vonMisesStress);

	const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[interval].material];
	bool rotateElasticity = mat->linear_elastic_properties.model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC;
	bool constCosSin = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
	// it is dirty hack just to be sure that compiler must assume both variants (currently settings.sigma = 0 and diffusion_split = false)
	bool constElasticity = !settings.contact_interfaces;
	bool constRotation = settings.load_steps == 1;
	if (mat->linear_elastic_properties.model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC) {
		if (mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN) {
			if (ndim == 2) {
				rotateElasticity &= mat->coordinate_system.rotation.z.isset;
			}
			if (ndim == 3) {
				rotateElasticity &= mat->coordinate_system.rotation.x.isset | mat->coordinate_system.rotation.y.isset | mat->coordinate_system.rotation.z.isset;
			}
		}
	}

	if (info::ecf->always_update_conductivity) { // TODO
		constElasticity = false;
	}
	bool storeCosSin = settings.reassembling_optimization && action == ActionOperator::ASSEMBLE;
	bool loadCosSin  = settings.reassembling_optimization && action != ActionOperator::ASSEMBLE;

	bool hasAcceleration = false;
	bool constAcceleration = true;
	auto AccelerationEval = configuration.acceleration.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (AccelerationEval != configuration.acceleration.end()) {
		hasAcceleration = true;
		constAcceleration = AccelerationEval->second.x.evaluator != nullptr && AccelerationEval->second.y.evaluator != nullptr && AccelerationEval->second.z.evaluator != nullptr;
	}

	bool hasVelocity = false;
	bool constVelocity = true;
	auto VelocityEval = configuration.angular_velocity.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (VelocityEval != configuration.angular_velocity.end()) {
		hasVelocity = true;
		constVelocity = VelocityEval->second.x.evaluator != nullptr && VelocityEval->second.y.evaluator != nullptr && VelocityEval->second.z.evaluator != nullptr;
	}

	bool constThickness = true;
	auto thicknessEval = settings.thickness.find(info::mesh->elementsRegions[info::mesh->elements->eintervals[interval].region]->name);
	if (thicknessEval != settings.thickness.end()) {
		constThickness = thicknessEval->second.evaluator != nullptr;
	}

	bool cooToGP = mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN || hasVelocity || axisymmetric;
	bool computeK = action == ActionOperator::ASSEMBLE || action == ActionOperator::REASSEMBLE;
	bool computeStress = action == ActionOperator::SOLUTION && info::ecf->output.results_selection.stress;
	bool computeElasticity = computeK | computeStress;
	bool getDisplacement = computeStress;

	esint chunks = elements / SIMD::size;
	double init = eslog::time() - initStart;

	if (cooToGP) {
		cooAndGps.move(SIMD::size);
	} else {
		coo.move(SIMD::size);
	}
	if (getDisplacement) {
		displacement.move(SIMD::size);
	}
	if (computeK) {
		stiffness.move(SIMD::size);
		if (hasAcceleration) {
			acceleration.move(SIMD::size);
		}
		if (hasVelocity) {
			velocity.move(SIMD::size);
		}
	}
	if (computeStress) {
		stress.move(SIMD::size);
	}

	double start = eslog::time();
	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xFACE); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xCAFE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xFEED); break; // TODO
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
		if (getDisplacement) {
			displacement.simd(element);
//			printf(" displacement");
		}
		if (!constThickness) {
			updateThickness(element, thicknessEval->second.evaluator);
//			printf(" thickness");
		}
		integration.simd(element);
//		printf(" integration");
//		if (getTemp) {
//			temp.simd(element);
//		}

		if (computeElasticity) {
			if (!constElasticity) {
				updateElasticity<DataDescriptor, nodes, gps, ndim, edim, etype>()(element, mat);
//				printf(" updateConductivity");
			}
			if (rotateElasticity) {
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
			stiffness.simd(element);
//			printf(" stiffness");
			if (hasAcceleration) {
				if (!constAcceleration) {
					updateAcc(element, AccelerationEval->second.x.evaluator);
//					printf(" updateAcc");
				}
				acceleration.simd(element);
//				printf(" acc");
			}
			if (hasVelocity) {
				if (!constVelocity) {
					updateVelocity(element, VelocityEval->second.x.evaluator);
//					printf(" updateVel");
				}
				velocity.simd(element);
//				printf(" velocity");
			}
		}
		if (computeStress) {
			stress.simd(element);
//			printf(" stress");
		}
//		printf("\n");
	}

	switch (action) {
	case ActionOperator::ASSEMBLE  : __SSC_MARK(0xDEAD); break;
	case ActionOperator::REASSEMBLE: __SSC_MARK(0xFADE); break;
	case ActionOperator::SOLUTION  : __SSC_MARK(0xBEED); break; // TODO
	default:
		eslog::error("unsupported action\n");
	}
	double loop = eslog::time() - start;

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
	return measurements(init, loop);
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateConditions2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return conditionsloop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return conditionsloop<StructuralMechanicsDataDescriptor, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return conditionsloop<StructuralMechanicsDataDescriptor, 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return conditionsloop<StructuralMechanicsDataDescriptor, 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype>(action, ops, interval, elements); break;
	default: return measurements();
	}
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiateConditions3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return conditionsloop<StructuralMechanicsDataDescriptor,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return conditionsloop<StructuralMechanicsDataDescriptor, 10, StructuralMechanicsGPC::TETRA10   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return conditionsloop<StructuralMechanicsDataDescriptor,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return conditionsloop<StructuralMechanicsDataDescriptor, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return conditionsloop<StructuralMechanicsDataDescriptor,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return conditionsloop<StructuralMechanicsDataDescriptor, 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return conditionsloop<StructuralMechanicsDataDescriptor,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, etype>(action, ops, interval, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return conditionsloop<StructuralMechanicsDataDescriptor, 20, StructuralMechanicsGPC::HEXA20    , 3, 3, etype>(action, ops, interval, elements); break;
	default: return measurements();
	}
}

Assembler::measurements StructuralMechanics::instantiateConditions(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return instantiateConditions2D<StructuralMechanicsElementType::SYMMETRIC_PLANE             >(action, code, ops, interval, elements);
		case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return instantiateConditions2D<StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_VOLUME: return instantiateConditions3D<StructuralMechanicsElementType::SYMMETRIC_VOLUME>(action, code, ops, interval, elements);
		}
	}
	return measurements();
}

}

