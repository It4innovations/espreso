
#include "structuralmechanics3d.elasticity.kernel.h"
#include "heattransfer3d.kernel.h"
#include "esinfo/stepinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/preprocessing/meshpreprocessing.h"
#include "physics/system/builder/builder.h"
#include "basis/containers/point.h"
#include "basis/evaluator/evaluator.h"
#include "basis/utilities/parser.h"
#include "config/ecf/physics/heattransfer.h"
#include "math/matrix.dense.h"
#include "math/vector.dense.h"

#include <cmath>

using namespace espreso;

StructuralMechanics3DKernel::StructuralMechanics3DKernel(StructuralMechanics3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanics3DBaseKernel(previous, physics, gsettings, configuration)
{
	solutions.push_back(VectorDense(iterator.displacement.output.data->data.size(), iterator.displacement.output.data->data.data()));
	orientation = NULL;
	for (size_t i = 0; i < info::mesh->elements->data.size(); ++i) {
		if (StringCompare::caseInsensitiveEq("ORIENTATION", info::mesh->elements->data[i]->name)) {
			orientation = info::mesh->elements->data[i];
		}
	}
}

StructuralMechanics3DKernel::StructuralMechanics3DKernel(HeatTransfer3DKernel *previous, PhysicsConfiguration &physics, StructuralMechanicsGlobalSettings &gsettings, StructuralMechanicsLoadStepConfiguration &configuration)
: StructuralMechanics3DBaseKernel(previous, physics, gsettings, configuration)
{
	solutions.push_back(VectorDense(iterator.displacement.output.data->data.size(), iterator.displacement.output.data->data.data()));
	orientation = NULL;
	for (size_t i = 0; i < info::mesh->elements->data.size(); ++i) {
		if (StringCompare::caseInsensitiveEq("ORIENTATION", info::mesh->elements->data[i]->name)) {
			orientation = info::mesh->elements->data[i];
		}
	}
}

StructuralMechanics3DKernel::~StructuralMechanics3DKernel()
{

}

void StructuralMechanics3DKernel::assembleLinearElasticMaterialMatrix(esint node, double *coordinates, const MaterialBaseConfiguration *mat, double time, double temp, MatrixDense &K, Point &orientation) const
{
	double Ex, Ey, Ez, miXY, miXZ, miYZ, Gx, Gy, Gz;
	Point p(coordinates[0], coordinates[1], coordinates[2]);

	Evaluator::Params params;
	params.coords(3, coordinates);
	params.temp(&temp);

	MatrixDense TCT(6, 6), T(6, 6), C(6, 6), CT(6, 6);

	switch (mat->linear_elastic_properties.model) {

	case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC: {
		Ex = Ey = Ez = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->eval(params);
		miXY = miXZ = miYZ = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->eval(params);

		double EE = Ex / ((1 + miXY) * (1 - 2 * miXY));

		C(0, 0) = EE * (1.0 - miXY);
		C(1, 1) = EE * (1.0 - miXY);
		C(2, 2) = EE * (1.0 - miXY);
		C(3, 3) = EE * (0.5 - miXY);
		C(4, 4) = EE * (0.5 - miXY);
		C(5, 5) = EE * (0.5 - miXY);

		C(0, 1) = EE * miXY;
		C(0, 2) = EE * miXY;
		C(0, 3) = 0;
		C(0, 4) = 0;
		C(0, 5) = 0;
		C(1, 2) = EE * miXY;
		C(1, 3) = 0;
		C(1, 4) = 0;
		C(1, 5) = 0;
		C(2, 3) = 0;
		C(2, 4) = 0;
		C(2, 5) = 0;
		C(3, 4) = 0;
		C(3, 5) = 0;
		C(4, 5) = 0;

		C(1, 0) = EE * miXY;
		C(2, 0) = EE * miXY;
		C(2, 1) = EE * miXY;
		C(3, 0) = 0;
		C(3, 1) = 0;
		C(3, 2) = 0;
		C(4, 0) = 0;
		C(4, 1) = 0;
		C(4, 2) = 0;
		C(4, 3) = 0;
		C(5, 0) = 0;
		C(5, 1) = 0;
		C(5, 2) = 0;
		C(5, 3) = 0;
		C(5, 4) = 0;
	} break;

	case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC: {
		C(0, 0) = mat->linear_elastic_properties.anisotropic.get(0, 0).evaluator->eval(params);
		C(1, 1) = mat->linear_elastic_properties.anisotropic.get(1, 1).evaluator->eval(params);
		C(2, 2) = mat->linear_elastic_properties.anisotropic.get(2, 2).evaluator->eval(params);
		C(3, 3) = mat->linear_elastic_properties.anisotropic.get(3, 3).evaluator->eval(params);
		C(4, 4) = mat->linear_elastic_properties.anisotropic.get(4, 4).evaluator->eval(params);
		C(5, 5) = mat->linear_elastic_properties.anisotropic.get(5, 5).evaluator->eval(params);

		C(0, 1) = mat->linear_elastic_properties.anisotropic.get(0, 1).evaluator->eval(params);
		C(0, 2) = mat->linear_elastic_properties.anisotropic.get(0, 2).evaluator->eval(params);
		C(0, 3) = mat->linear_elastic_properties.anisotropic.get(0, 3).evaluator->eval(params);
		C(0, 4) = mat->linear_elastic_properties.anisotropic.get(0, 4).evaluator->eval(params);
		C(0, 5) = mat->linear_elastic_properties.anisotropic.get(0, 5).evaluator->eval(params);
		C(1, 2) = mat->linear_elastic_properties.anisotropic.get(1, 2).evaluator->eval(params);
		C(1, 3) = mat->linear_elastic_properties.anisotropic.get(1, 3).evaluator->eval(params);
		C(1, 4) = mat->linear_elastic_properties.anisotropic.get(1, 4).evaluator->eval(params);
		C(1, 5) = mat->linear_elastic_properties.anisotropic.get(1, 5).evaluator->eval(params);
		C(2, 3) = mat->linear_elastic_properties.anisotropic.get(2, 3).evaluator->eval(params);
		C(2, 4) = mat->linear_elastic_properties.anisotropic.get(2, 4).evaluator->eval(params);
		C(2, 5) = mat->linear_elastic_properties.anisotropic.get(2, 5).evaluator->eval(params);
		C(3, 4) = mat->linear_elastic_properties.anisotropic.get(3, 4).evaluator->eval(params);
		C(3, 5) = mat->linear_elastic_properties.anisotropic.get(3, 5).evaluator->eval(params);
		C(4, 5) = mat->linear_elastic_properties.anisotropic.get(4, 5).evaluator->eval(params);

		C(1, 0) = C(0, 1);
		C(2, 0) = C(0, 2);
		C(2, 1) = C(1, 2);
		C(3, 0) = C(0, 3);
		C(3, 1) = C(1, 3);
		C(3, 2) = C(2, 3);
		C(4, 0) = C(0, 4);
		C(4, 1) = C(1, 4);
		C(4, 2) = C(2, 4);
		C(4, 3) = C(3, 4);
		C(5, 0) = C(0, 5);
		C(5, 1) = C(1, 5);
		C(5, 2) = C(2, 5);
		C(5, 3) = C(3, 5);
		C(5, 4) = C(4, 5);
	} break;

	case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC: {
		Ex = mat->linear_elastic_properties.young_modulus.get(0, 0).evaluator->eval(params);
		Ey = mat->linear_elastic_properties.young_modulus.get(1, 1).evaluator->eval(params);
		Ez = mat->linear_elastic_properties.young_modulus.get(2, 2).evaluator->eval(params);

		miXY = mat->linear_elastic_properties.poisson_ratio.get(0, 0).evaluator->eval(params);
		miXZ = mat->linear_elastic_properties.poisson_ratio.get(1, 1).evaluator->eval(params);
		miYZ = mat->linear_elastic_properties.poisson_ratio.get(2, 2).evaluator->eval(params);

		Gx = mat->linear_elastic_properties.shear_modulus.get(0, 0).evaluator->eval(params);
		Gy = mat->linear_elastic_properties.shear_modulus.get(1, 1).evaluator->eval(params);
		Gz = mat->linear_elastic_properties.shear_modulus.get(2, 2).evaluator->eval(params);

		double miYX = miXY * Ey / Ex;
		double miZY = miYZ * Ez / Ey;
		double miZX = miXZ * Ex / Ez;

		double ksi = 1 - (miXY * miYX + miYZ * miZY + miXZ * miZX) - (miXY * miYZ * miZX + miYX * miZY * miXZ);

		double dxx = Ex * (1 - miYZ * miZY) / ksi;
		double dxy = Ey * (miXY + miXZ * miZY) / ksi;
		double dxz = Ez * (miXZ + miYZ * miXY)  /ksi;
		double dyy = Ey * (1 - miXZ * miZX) / ksi;
		double dyz = Ez * (miYZ + miYX * miXZ) / ksi;
		double dzz = Ez * (1 - miYX * miXY) / ksi;

		C(0, 0) = dxx;
		C(1, 1) = dyy;
		C(2, 2) = dzz;
		C(3, 3) = Gx;
		C(4, 4) = Gz;
		C(5, 5) = Gy;

		C(0, 1) = dxy;
		C(0, 2) = dxz;
		C(0, 3) = 0;
		C(0, 4) = 0;
		C(0, 5) = 0;
		C(1, 2) = dyz;
		C(1, 3) = 0;
		C(1, 4) = 0;
		C(1, 5) = 0;
		C(2, 3) = 0;
		C(2, 4) = 0;
		C(2, 5) = 0;
		C(3, 4) = 0;
		C(3, 5) = 0;
		C(4, 5) = 0;

		C(1, 0) = dxy;
		C(2, 0) = dxz;
		C(2, 1) = dyz;
		C(3, 0) = 0;
		C(3, 1) = 0;
		C(3, 2) = 0;
		C(4, 0) = 0;
		C(4, 1) = 0;
		C(4, 2) = 0;
		C(4, 3) = 0;
		C(5, 0) = 0;
		C(5, 1) = 0;
		C(5, 2) = 0;
		C(5, 3) = 0;
		C(5, 4) = 0;
	} break;

	default:
		eslog::error("Structural mechanics 3D not supports set material model.\n");
	}

	if (mat->linear_elastic_properties.orientation) {
		Point angle = -orientation * (M_PI / 180);

		T[0][0] = std::pow(std::cos(angle.x), 2.0) * std::pow(std::cos(angle.z), 2.0) - std::sin(angle.x * 2.0) * std::sin(angle.z * 2.0) * std::cos(angle.y) * (1.0 / 2.0) + std::pow(std::cos(angle.y), 2.0) * std::pow(std::sin(angle.x), 2.0) * std::pow(std::sin(angle.z), 2.0);
		T[0][1] = std::pow(std::cos(angle.z), 2.0) * std::pow(std::sin(angle.x), 2.0) + std::sin(angle.x * 2.0) * std::sin(angle.z * 2.0) * std::cos(angle.y) * (1.0 / 2.0) + std::pow(std::cos(angle.x), 2.0) * std::pow(std::cos(angle.y), 2.0) * std::pow(std::sin(angle.z), 2.0);
		T[0][2] = std::pow(std::sin(angle.y), 2.0) * std::pow(std::sin(angle.z), 2.0);
		T[0][3] = std::sin(angle.z * 2.0) * std::sin(angle.x) * std::sin(angle.y) + std::sin(angle.y * 2.0) * std::cos(angle.x) * std::pow(std::sin(angle.z), 2.0);
		T[0][4] = std::sin(angle.z * 2.0) * std::cos(angle.x) * std::sin(angle.y) - std::sin(angle.y * 2.0) * std::sin(angle.x) * std::pow(std::sin(angle.z), 2.0);
		T[0][5] = std::sin(angle.x * 2.0) * std::pow(std::cos(angle.z), 2.0) + std::cos(angle.x * 2.0) * std::sin(angle.z * 2.0) * std::cos(angle.y) - std::sin(angle.x * 2.0) * std::pow(std::cos(angle.y), 2.0) * std::pow(std::sin(angle.z), 2.0);
		T[1][0] = std::pow(std::cos(angle.x), 2.0) * std::pow(std::sin(angle.z), 2.0) + std::sin(angle.x * 2.0) * std::sin(angle.z * 2.0) * std::cos(angle.y) * (1.0 / 2.0) + std::pow(std::cos(angle.y), 2.0) * std::pow(std::cos(angle.z), 2.0) * std::pow(std::sin(angle.x), 2.0);
		T[1][1] = std::pow(std::sin(angle.x), 2.0) * std::pow(std::sin(angle.z), 2.0) - std::sin(angle.x * 2.0) * std::sin(angle.z * 2.0) * std::cos(angle.y) * (1.0 / 2.0) + std::pow(std::cos(angle.x), 2.0) * std::pow(std::cos(angle.y), 2.0) * std::pow(std::cos(angle.z), 2.0);
		T[1][2] = std::pow(std::cos(angle.z), 2.0) * std::pow(std::sin(angle.y), 2.0);
		T[1][3] = -std::sin(angle.z * 2.0) * std::sin(angle.x) * std::sin(angle.y) + std::sin(angle.y * 2.0) * std::cos(angle.x) * std::pow(std::cos(angle.z), 2.0);
		T[1][4] = -std::sin(angle.z * 2.0) * std::cos(angle.x) * std::sin(angle.y) - std::sin(angle.y * 2.0) * std::pow(std::cos(angle.z), 2.0) * std::sin(angle.x);
		T[1][5] = std::sin(angle.x * 2.0) * std::pow(std::sin(angle.z), 2.0) - std::cos(angle.x * 2.0) * std::sin(angle.z * 2.0) * std::cos(angle.y) - std::sin(angle.x * 2.0) * std::pow(std::cos(angle.y), 2.0) * std::pow(std::cos(angle.z), 2.0);
		T[2][0] = std::pow(std::sin(angle.x), 2.0) * std::pow(std::sin(angle.y), 2.0);
		T[2][1] = std::pow(std::cos(angle.x), 2.0) * std::pow(std::sin(angle.y), 2.0);
		T[2][2] = std::pow(std::cos(angle.y), 2.0);
		T[2][3] = -std::sin(angle.y * 2.0) * std::cos(angle.x);
		T[2][4] = std::sin(angle.y * 2.0) * std::sin(angle.x);
		T[2][5] = -std::sin(angle.x * 2.0) * std::pow(std::sin(angle.y), 2.0);
		T[3][0] = std::sin(angle.x * 2.0) * std::sin(angle.y) * std::sin(angle.z) * (-1.0 / 2.0) - std::sin(angle.y * 2.0) * std::cos(angle.z) * std::pow(std::sin(angle.x), 2.0) * (1.0 / 2.0);
		T[3][1] = std::sin(angle.x * 2.0) * std::sin(angle.y) * std::sin(angle.z) * (1.0 / 2.0) - std::sin(angle.y * 2.0) * std::pow(std::cos(angle.x), 2.0) * std::cos(angle.z) * (1.0 / 2.0);
		T[3][2] = std::sin(angle.y * 2.0) * std::cos(angle.z) * (1.0 / 2.0);
		T[3][3] = std::cos(angle.y * 2.0) * std::cos(angle.x) * std::cos(angle.z) - std::cos(angle.y) * std::sin(angle.x) * std::sin(angle.z);
		T[3][4] = -std::cos(angle.y * 2.0) * std::cos(angle.z) * std::sin(angle.x) - std::cos(angle.x) * std::cos(angle.y) * std::sin(angle.z);
		T[3][5] = std::cos(angle.x * 2.0) * std::sin(angle.y) * std::sin(angle.z) + std::sin(angle.x * 2.0) * std::sin(angle.y * 2.0) * std::cos(angle.z) * (1.0 / 2.0);
		T[4][0] = std::sin(angle.x * 2.0) * std::cos(angle.z) * std::sin(angle.y) * (1.0 / 2.0) - std::sin(angle.y * 2.0) * std::pow(std::sin(angle.x), 2.0) * std::sin(angle.z) * (1.0 / 2.0);
		T[4][1] = std::sin(angle.x * 2.0) * std::cos(angle.z) * std::sin(angle.y) * (-1.0 / 2.0) - std::sin(angle.y * 2.0) * std::pow(std::cos(angle.x), 2.0) * std::sin(angle.z) * (1.0 / 2.0);
		T[4][2] = std::sin(angle.y * 2.0) * std::sin(angle.z) * (1.0 / 2.0);
		T[4][3] = std::cos(angle.y * 2.0) * std::cos(angle.x) * std::sin(angle.z) + std::cos(angle.y) * std::cos(angle.z) * std::sin(angle.x);
		T[4][4] = -std::cos(angle.y * 2.0) * std::sin(angle.x) * std::sin(angle.z) + std::cos(angle.x) * std::cos(angle.y) * std::cos(angle.z);
		T[4][5] = -std::cos(angle.x * 2.0) * std::cos(angle.z) * std::sin(angle.y) + std::sin(angle.x * 2.0) * std::sin(angle.y * 2.0) * std::sin(angle.z) * (1.0 / 2.0);
		T[5][0] = std::sin(angle.z * 2.0) * std::pow(std::cos(angle.x), 2.0) * (-1.0 / 2.0) - std::cos(angle.z * 2.0) * std::sin(angle.x * 2.0) * std::cos(angle.y) * (1.0 / 2.0) + std::sin(angle.z * 2.0) * std::pow(std::cos(angle.y), 2.0) * std::pow(std::sin(angle.x), 2.0) * (1.0 / 2.0);
		T[5][1] = std::sin(angle.z * 2.0) * std::pow(std::sin(angle.x), 2.0) * (-1.0 / 2.0) + std::cos(angle.z * 2.0) * std::sin(angle.x * 2.0) * std::cos(angle.y) * (1.0 / 2.0) + std::sin(angle.z * 2.0) * std::pow(std::cos(angle.x), 2.0) * std::pow(std::cos(angle.y), 2.0) * (1.0 / 2.0);
		T[5][2] = std::sin(angle.z * 2.0) * std::pow(std::sin(angle.y), 2.0) * (1.0 / 2.0);
		T[5][3] = std::cos(angle.z * 2.0) * std::sin(angle.x) * std::sin(angle.y) + std::sin(angle.y * 2.0) * std::sin(angle.z * 2.0) * std::cos(angle.x) * (1.0 / 2.0);
		T[5][4] = std::cos(angle.z * 2.0) * std::cos(angle.x) * std::sin(angle.y) - std::sin(angle.y * 2.0) * std::sin(angle.z * 2.0) * std::sin(angle.x) * (1.0 / 2.0);
		T[5][5] = std::sin(angle.x * 2.0) * std::sin(angle.z * 2.0) * (-1.0 / 2.0) + std::cos(angle.x * 2.0) * std::cos(angle.z * 2.0) * std::cos(angle.y) - std::sin(angle.x * 2.0) * std::sin(angle.z * 2.0) * std::pow(std::cos(angle.y), 2.0) * (1.0 / 2.0);

		CT.multiply(C, T, 1, 0, false, true);
		TCT.multiply(T, CT);
	} else {
		TCT = C;
	}

	size_t k = 0;
	for (size_t i = 0; i < 6; i++) {
		K(node, k++) = TCT(i, i);
	}
	for (size_t i = 0; i < 6; i++) {
		for (size_t j = i + 1; j < 6; j++) {
			K(node, k++) = TCT(i, j);
		}
	}
	for (size_t i = 0; i < 6; i++) {
		for (size_t j = 0; j < i; j++) {
			K(node, k++) = TCT(i, j);
		}
	}
}

void StructuralMechanics3DKernel::assembleHyperElasticMaterialMatrix(const MaterialBaseConfiguration *mat, MatrixDense &F, MatrixDense &C, MatrixDense &S) const
{
	C.fill(0);
	MatrixDense rCG(3, 3);
	if (!step::isInitial()) {
		rCG.multiply(F, F, 1, 0, true, false);
	}
	double c01, c10, d, dd, K, G, EE = 0, miXY = 0, lambdaL, E = 0;
	MatrixDense eHat(3, 3), eVec(6, 1), sVec(6, 1);

	if (step::isInitial()) {
		switch (mat->hyper_elastic_properties.model) {
		case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_CMP:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_INC:
			E = mat->hyper_elastic_properties.E.evaluator->eval(Evaluator::Params());
			miXY = mat->hyper_elastic_properties.mi.evaluator->eval(Evaluator::Params());
			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_2PARAMS:
			c01 = mat->hyper_elastic_properties.C01.evaluator->eval(Evaluator::Params());
			c10 = mat->hyper_elastic_properties.C10.evaluator->eval(Evaluator::Params());
			d = mat->hyper_elastic_properties.d.evaluator->eval(Evaluator::Params());

			G = 2 * (c10 + c01);
			dd = (1 - 2 * G) / (c10 + c01);
			K = 2 / dd;
			E = (9 * K * G) / (3 * K + G);
			miXY = (3 * K - 2 * G) / (2 * (3 * K + G));

			E = 2e11;
			miXY = 0.3;

			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_3PARAMS:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_5PARAMS:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_9PARAMS:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::ARRUDA_BOYCE:
			d = mat->hyper_elastic_properties.d.evaluator->eval(Evaluator::Params());
			G = mat->hyper_elastic_properties.G.evaluator->eval(Evaluator::Params());
			K = 2 / d;
			E  = (9 * K * G) / (3 * K + G);
			miXY = (3 * K - 2 * G) / (2 * (3 * K + G));
			break;
		case HyperElasticPropertiesConfiguration::MODEL::BLATZ_KO_FOAM:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::GENT:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::OGDEN_1:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::OGDEN_2:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::OGDEN_3:
			break;
		}

		EE = E / ((1 + miXY) * (1 - 2 * miXY));
		C(0, 0) = EE * (1.0 - miXY);
		C(1, 1) = EE * (1.0 - miXY);
		C(2, 2) = EE * (1.0 - miXY);
		C(3, 3) = EE * (0.5 - miXY);
		C(4, 4) = EE * (0.5 - miXY);
		C(5, 5) = EE * (0.5 - miXY);
		C(0, 1) = EE * miXY;
		C(0, 2) = EE * miXY;
		C(1, 2) = EE * miXY;
		C(1, 0) = EE * miXY;
		C(2, 0) = EE * miXY;
		C(2, 1) = EE * miXY;
	} else {
		double detF, detCG, lnDetF, lambda, mu;
		MatrixDense invCG(3, 3);

		switch (mat->hyper_elastic_properties.model) {
		case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_CMP:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::NEO_HOOKEN_INC:
			E = mat->hyper_elastic_properties.E.evaluator->eval(Evaluator::Params());
			miXY = mat->hyper_elastic_properties.mi.evaluator->eval(Evaluator::Params());
			lambda = E * miXY / ((1 + miXY) * (1 - 2 * miXY));
			mu = E / (2 * (1 + miXY));

			S(0,0) = mu+(mu*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2)))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2)))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2));
			S(0,1) = (mu*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0)-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0);
			S(0,2) = (mu*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0)-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0);
			S(1,0) = (mu*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0)-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0);
			S(1,1) = mu+(mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2)))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2)))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2));
			S(1,2) = (mu*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0)-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0);
			S(2,0) = (mu*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0)-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0);
			S(2,1) = (mu*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0)-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))*2.0+(rCG(0,2)*rCG(0,2))*rCG(1,1)*2.0+(rCG(0,1)*rCG(0,1))*rCG(2,2)*2.0-rCG(0,1)*rCG(0,2)*rCG(1,2)*4.0-rCG(0,0)*rCG(1,1)*rCG(2,2)*2.0);
			S(2,2) = mu+(mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1)))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1)))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2));

			C(0,0) = lambda*pow(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*pow(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*pow(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(0,1) = (mu*rCG(2,2)*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-(lambda*rCG(2,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(0,2) = (mu*rCG(1,1)*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-(lambda*rCG(1,1)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(0,3) = (lambda*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(0,4) = (mu*rCG(1,2)*-2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+(lambda*rCG(1,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(0,5) = (lambda*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(1,0) = (mu*rCG(2,2)*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-(lambda*rCG(2,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(1,1) = lambda*pow(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*pow(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*pow(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(1,2) = (mu*rCG(0,0)*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-(lambda*rCG(0,0)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(1,3) = (lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(1,4) = (lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(1,5) = (mu*rCG(0,2)*-2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+(lambda*rCG(0,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(2,0) = (mu*rCG(1,1)*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-(lambda*rCG(1,1)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(2,1) = (mu*rCG(0,0)*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-(lambda*rCG(0,0)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(2,2) = lambda*pow(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+mu*pow(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*pow(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1),2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)*2.0;
			C(2,3) = (mu*rCG(0,1)*-2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+(lambda*rCG(0,1)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(2,4) = (lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(2,5) = (lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(3,0) = (lambda*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(3,1) = (lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(3,2) = (mu*rCG(0,1)*-2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+(lambda*rCG(0,1)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(3,3) = (lambda*pow(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*pow(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(mu*rCG(2,2))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*pow(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+(lambda*rCG(2,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2));
			C(3,4) = (mu*rCG(0,2))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(lambda*rCG(0,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0;
			C(3,5) = (mu*rCG(1,2))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(lambda*rCG(1,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0;
			C(4,0) = (mu*rCG(1,2)*-2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+(lambda*rCG(1,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(4,1) = (lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(4,2) = (lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(4,3) = (mu*rCG(0,2))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(lambda*rCG(0,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0;
			C(4,4) = (lambda*pow(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*pow(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(mu*rCG(0,0))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*pow(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+(lambda*rCG(0,0)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2));
			C(4,5) = (mu*rCG(0,1))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(lambda*rCG(0,1)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0;
			C(5,0) = (lambda*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(1,1)*rCG(2,2)-rCG(1,2)*rCG(1,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(5,1) = (mu*rCG(0,2)*-2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)+(lambda*rCG(0,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*2.0)/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(2,2)-rCG(0,2)*rCG(0,2))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(5,2) = (lambda*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+mu*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0)-lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,0)*rCG(1,1)-rCG(0,1)*rCG(0,1))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0);
			C(5,3) = (mu*rCG(1,2))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(lambda*rCG(1,2)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*(rCG(0,2)*rCG(1,2)*2.0-rCG(0,1)*rCG(2,2)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0;
			C(5,4) = (mu*rCG(0,1))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))+(lambda*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(lambda*rCG(0,1)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*(rCG(0,1)*rCG(0,2)*2.0-rCG(0,0)*rCG(1,2)*2.0)*(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0;
			C(5,5) = (lambda*pow(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/4.0+(mu*pow(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0-(mu*rCG(1,1))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2))-(lambda*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2)))*pow(rCG(0,1)*rCG(1,2)*2.0-rCG(0,2)*rCG(1,1)*2.0,2.0)*1.0/pow(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2),2.0))/2.0+(lambda*rCG(1,1)*log(sqrt(-rCG(0,0)*(rCG(1,2)*rCG(1,2))-(rCG(0,2)*rCG(0,2))*rCG(1,1)-(rCG(0,1)*rCG(0,1))*rCG(2,2)+rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0+rCG(0,0)*rCG(1,1)*rCG(2,2))))/(rCG(0,0)*(rCG(1,2)*rCG(1,2))+(rCG(0,2)*rCG(0,2))*rCG(1,1)+(rCG(0,1)*rCG(0,1))*rCG(2,2)-rCG(0,1)*rCG(0,2)*rCG(1,2)*2.0-rCG(0,0)*rCG(1,1)*rCG(2,2));
			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_2PARAMS:
			c01 = mat->hyper_elastic_properties.C01.evaluator->eval(Evaluator::Params());
			c10 = mat->hyper_elastic_properties.C10.evaluator->eval(Evaluator::Params());
			d = mat->hyper_elastic_properties.d.evaluator->eval(Evaluator::Params());

			K = 2 / d;
			G = 2 * (c10 + c01);

			detF = MATH::determinant3x3(F.vals);
			detCG = MATH::determinant3x3(rCG.vals);
			MATH::Dense3x3inverse(F.vals, invCG.vals, detCG);

			lnDetF = log(detF);

			lambda = 17000;
			mu = 1e5;
			S(0, 0) = lambda * lnDetF * invCG(0, 0) + mu * (1 - invCG(0, 0));
			S(1, 1) = lambda * lnDetF * invCG(1, 1) + mu * (1 - invCG(1, 1));
			S(2, 2) = lambda * lnDetF * invCG(2, 2) + mu * (1 - invCG(2, 2));

			S(0, 1) = lambda * lnDetF * invCG(0, 1) + mu * (0 - invCG(0, 1));
			S(0, 2) = lambda * lnDetF * invCG(0, 2) + mu * (0 - invCG(0, 2));
			S(1, 0) = lambda * lnDetF * invCG(1, 0) + mu * (0 - invCG(1, 0));
			S(1, 2) = lambda * lnDetF * invCG(1, 2) + mu * (0 - invCG(1, 2));
			S(2, 0) = lambda * lnDetF * invCG(2, 0) + mu * (0 - invCG(2, 0));
			S(2, 1) = lambda * lnDetF * invCG(2, 1) + mu * (0 - invCG(2, 1));

			C(1 - 1, 1 - 1) = (-2 * lnDetF * lambda+lambda+2*mu)*invCG(0, 0) * invCG(0, 0);
			C(1 - 1, 2 - 1) = (2*mu-2*lnDetF*lambda)*invCG(0, 1)*invCG(0, 1) + lambda*invCG(0, 0)*invCG(1, 1);
			C(1 - 1, 3 - 1) = (2*mu-2*lnDetF*lambda)*invCG(2, 0)*invCG(2, 0) + lambda*invCG(0, 0)*invCG(2, 2);
			C(1 - 1, 4 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(0, 1)*invCG(0, 0);
			C(1 - 1, 5 - 1) = lambda*invCG(0, 0)*invCG(1, 2)+2*(mu-lnDetF*lambda)*invCG(0, 1)*invCG(2, 0);
			C(1 - 1, 6 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(0, 0)*invCG(2, 0);
			C(2 - 1, 2 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(1, 1)*invCG(1, 1);
			C(2 - 1, 3 - 1) = (2*mu-2*lnDetF*lambda)*invCG(1, 2)*invCG(1, 2) + lambda*invCG(2, 2)*invCG(1, 1);
			C(2 - 1, 4 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(0, 1)*invCG(1, 1);
			C(2 - 1, 5 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(1, 1)*invCG(1, 2);
			C(2 - 1, 6 - 1) = 2*(mu-lnDetF*lambda)*invCG(0, 1)*invCG(1, 2)+lambda*invCG(1, 1)*invCG(2, 0);
			C(3 - 1, 3 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(2, 2)*invCG(2, 2);
			C(3 - 1, 4 - 1) = 2*(mu-lnDetF*lambda)*invCG(1, 2)*invCG(2, 0)+lambda*invCG(0, 1)*invCG(2, 2);
			C(3 - 1, 5 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(1, 2)*invCG(2, 2);
			C(3 - 1, 6 - 1) = (-2*lnDetF*lambda+lambda+2*mu)*invCG(2, 0)*invCG(2, 2);
			C(4 - 1, 4 - 1) = (-lnDetF*lambda+lambda+mu)*invCG(0, 1)*invCG(0, 1)+(mu-lnDetF*lambda)*invCG(0, 0)*invCG(1, 1);
			C(4 - 1, 5 - 1) = (-lnDetF*lambda+lambda+mu)*invCG(0, 1)*invCG(1, 2)+(mu-lnDetF*lambda)*invCG(1, 1)*invCG(2, 0);
			C(4 - 1, 6 - 1) = (-lnDetF*lambda+lambda+mu)*invCG(0, 1)*invCG(2, 0)+(mu-lnDetF*lambda)*invCG(0, 0)*invCG(1, 2);
			C(5 - 1, 5 - 1) = (-lnDetF*lambda+lambda+mu)*invCG(1, 2)*invCG(1, 2)+(mu-lnDetF*lambda)*invCG(1, 1)*invCG(2, 2);
			C(5 - 1, 6 - 1) = (-lnDetF*lambda+lambda+mu)*invCG(1, 2)*invCG(2, 0)+(mu-lnDetF*lambda)*invCG(0, 1)*invCG(2, 2);
			C(6 - 1, 6 - 1) = (-lnDetF*lambda+lambda+mu)*invCG(2, 0)*invCG(2, 0)+(mu-lnDetF*lambda)*invCG(0, 0)*invCG(2, 2);


//			E = (9 * K * G) / (3 * K + G);
//			miXY = (3 * K - 2 * G) / (2 * (3 * K + G));

//			E = 2e11;
//			miXY = 0.3;
//
//			EE = E / ((1 + miXY) * (1 - 2 * miXY));
//			C(0, 0) = EE * (1.0 - miXY);
//			C(1, 1) = EE * (1.0 - miXY);
//			C(2, 2) = EE * (1.0 - miXY);
//			C(3, 3) = EE * (0.5 - miXY);
//			C(4, 4) = EE * (0.5 - miXY);
//			C(5, 5) = EE * (0.5 - miXY);
//			C(0, 1) = EE * miXY;
//			C(0, 2) = EE * miXY;
//			C(1, 2) = EE * miXY;
//			C(1, 0) = EE * miXY;
//			C(2, 0) = EE * miXY;
//			C(2, 1) = EE * miXY;
//
//			eHat.multiply(F, F, .5, 0, true, false);
//			eHat(0, 0) -= .5;
//			eHat(1, 1) -= .5;
//			eHat(2, 2) -= .5;
//			eVec(0, 0) = eHat(0, 0);
//			eVec(1, 0) = eHat(1, 1);
//			eVec(2, 0) = eHat(2, 2);
//			eVec(3, 0) = 2 * eHat(0, 1);
//			eVec(4, 0) = 2 * eHat(1, 2);
//			eVec(5, 0) = 2 * eHat(0, 2);
//			sVec.multiply(C, eVec);
//			S(0, 0) = sVec(0, 0);
//			S(1, 1) = sVec(1, 0);
//			S(2, 2) = sVec(2, 0);
//			S(0, 1) = sVec(3, 0);
//			S(1, 0) = sVec(3, 0);
//			S(0, 2) = sVec(5, 0);
//			S(2, 0) = sVec(5, 0);
//			S(1, 2) = sVec(4, 0);
//			S(2, 1) = sVec(4, 0);
//
//			S(0, 0) = c10+c01*(rCG(1, 1)+rCG(2, 2))+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/d;
//			S(0, 1) = -c01*rCG(0, 1)+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/(d*2.0);
//			S(0, 2) = -c01*rCG(0, 2)+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/(d*2.0);
//			S(1, 0) = -c01*rCG(0, 1)+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/(d*2.0);
//			S(1, 1) = c10+c01*(rCG(0, 0)+rCG(2, 2))+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/d;
//			S(1, 2) = -c01*rCG(1, 2)+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/(d*2.0);
//			S(2, 0) = -c01*rCG(0, 2)+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/(d*2.0);
//			S(2, 1) = -c01*rCG(1, 2)+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/(d*2.0);
//			S(2, 2) = c10+c01*(rCG(0, 0)+rCG(1, 1))+((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/d;
//
//			C(0, 0) = (pow(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2),2.0)*-2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*pow(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2),2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(0, 1) = c01*4.0+(rCG(2, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*4.0)/d-((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(0, 2) = c01*4.0+(rCG(1, 1)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*4.0)/d-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(0, 3) = -((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(0, 4) = (rCG(1, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*-4.0)/d-((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(0, 5) = -((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(1, 0) = c01*4.0+(rCG(2, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*4.0)/d-((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(1, 1) = (pow(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2),2.0)*-2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*pow(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2),2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(1, 2) = c01*4.0+(rCG(0, 0)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*4.0)/d-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(1, 3) = -((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(1, 4) = -((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(1, 5) = (rCG(0, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*-4.0)/d-((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(2, 0) = c01*4.0+(rCG(1, 1)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*4.0)/d-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(2, 1) = c01*4.0+(rCG(0, 0)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*4.0)/d-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(2, 2) = (pow(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1),2.0)*-2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*pow(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1),2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0)*2.0)/d;
//			C(2, 3) = (rCG(0, 1)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*-4.0)/d-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(2, 4) = -((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(2, 5) = -((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(3, 0) = -((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(3, 1) = -((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(3, 2) = (rCG(0, 1)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*-4.0)/d-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(3, 3) = c01*-2.0-pow(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0,2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)-(rCG(2, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*pow(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0,2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(3, 4) = ((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*(-1.0/2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))+(rCG(0, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(3, 5) = ((rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*(-1.0/2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))+(rCG(1, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(4, 0) = (rCG(1, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*-4.0)/d-((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(4, 1) = -((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(4, 2) = -((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(4, 3) = ((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*(-1.0/2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))+(rCG(0, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(4, 4) = c01*-2.0-pow(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0,2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)-(rCG(0, 0)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*pow(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0,2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(4, 5) = ((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(-1.0/2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))+(rCG(0, 1)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(5, 0) = -((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(5, 1) = (rCG(0, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*-4.0)/d-((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(5, 2) = -((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/d;
//			C(5, 3) = ((rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*(-1.0/2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))+(rCG(1, 2)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(5, 4) = ((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(-1.0/2.0))/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))+(rCG(0, 1)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);
//			C(5, 5) = c01*-2.0-pow(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0,2.0)/(d*(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)-(rCG(1, 1)*(sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*1.0/sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))*2.0)/d-((sqrt(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-1.0)*pow(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0,2.0)*1.0/pow(-rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))-(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)-(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)+rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0+rCG(0, 0)*rCG(1, 1)*rCG(2, 2),3.0/2.0))/(d*2.0);

			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_3PARAMS:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_5PARAMS:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::MOONEY_RIVLIN_9PARAMS:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::ARRUDA_BOYCE:
			d = mat->hyper_elastic_properties.d.evaluator->eval(Evaluator::Params());
			G = mat->hyper_elastic_properties.G.evaluator->eval(Evaluator::Params());
			lambdaL = mat->hyper_elastic_properties.lambdaL.evaluator->eval(Evaluator::Params());

			S(0, 0) = ((rCG(1, 1)*rCG(2, 2))/2.0+(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)-(rCG(1, 2)*rCG(1, 2))/2.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*(1.1E1/3.5E2)+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.078014184397163E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),4.0)*3.851576994434137E-3+(1.0/(lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0))/2.0E1+1.0/2.0);
			S(0, 1) = (rCG(0, 2)*rCG(1, 2)-rCG(0, 1)*rCG(2, 2)+(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0))/(d*2.0);
			S(0, 2) = (rCG(0, 1)*rCG(1, 2)-rCG(0, 2)*rCG(1, 1)+(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0))/(d*2.0);
			S(1, 0) = (rCG(0, 2)*rCG(1, 2)-rCG(0, 1)*rCG(2, 2)+(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0))/(d*2.0);
			S(1, 1) = ((rCG(0, 0)*rCG(2, 2))/2.0+(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)-(rCG(0, 2)*rCG(0, 2))/2.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*(1.1E1/3.5E2)+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.078014184397163E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),4.0)*3.851576994434137E-3+(1.0/(lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0))/2.0E1+1.0/2.0);
			S(1, 2) = (rCG(0, 1)*rCG(0, 2)-rCG(0, 0)*rCG(1, 2)+(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0))/(d*2.0);
			S(2, 0) = (rCG(0, 1)*rCG(1, 2)-rCG(0, 2)*rCG(1, 1)+(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0))/(d*2.0);
			S(2, 1) = (rCG(0, 1)*rCG(0, 2)-rCG(0, 0)*rCG(1, 2)+(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0))/(d*2.0);
			S(2, 2) = ((rCG(0, 0)*rCG(1, 1))/2.0+(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)-(rCG(0, 1)*rCG(0, 1))/2.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*(1.1E1/3.5E2)+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.078014184397163E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),4.0)*3.851576994434137E-3+(1.0/(lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0))/2.0E1+1.0/2.0);

			C(0, 0) = G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0+(pow(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2),2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0)*2.0)/d;
			C(0, 1) = ((rCG(2, 2)/2.0+rCG(2, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)+((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*4.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0;
			C(0, 2) = ((rCG(1, 1)/2.0+rCG(1, 1)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)+((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*4.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0;
			C(0, 3) = ((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(0, 4) = ((rCG(1, 2)+rCG(1, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*-2.0)/d;
			C(0, 5) = ((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(1, 0) = ((rCG(2, 2)/2.0+rCG(2, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)+((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*4.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0;
			C(1, 1) = G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0+(pow(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2),2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0)*2.0)/d;
			C(1, 2) = ((rCG(0, 0)/2.0+rCG(0, 0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)+((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*4.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0;
			C(1, 3) = ((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(1, 4) = ((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(1, 5) = ((rCG(0, 2)+rCG(0, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*-2.0)/d;
			C(2, 0) = ((rCG(1, 1)/2.0+rCG(1, 1)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)+((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*4.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0;
			C(2, 1) = ((rCG(0, 0)/2.0+rCG(0, 0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))*2.0+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)*2.0+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)*2.0-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*4.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)*2.0)+((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*4.0)/d+G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0;
			C(2, 2) = G*(1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),2.0)*3.234042553191489E-2+1.0/(lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL*lambdaL)*pow(rCG(0, 0)+rCG(1, 1)+rCG(2, 2),3.0)*1.540630797773655E-2+1.0/(lambdaL*lambdaL)/1.0E1+1.0/(lambdaL*lambdaL*lambdaL*lambdaL)*(rCG(0, 0)*2.0+rCG(1, 1)*2.0+rCG(2, 2)*2.0)*(1.1E1/3.5E2))*4.0+(pow(rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1),2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0)*2.0)/d;
			C(2, 3) = ((rCG(0, 1)+rCG(0, 1)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*-2.0)/d;
			C(2, 4) = ((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(2, 5) = ((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(3, 0) = ((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(3, 1) = ((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(3, 2) = ((rCG(0, 1)+rCG(0, 1)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*-2.0)/d;
			C(3, 3) = -(rCG(2, 2)-(pow(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0,2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0+rCG(2, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/d;
			C(3, 4) = (rCG(0, 2)+rCG(0, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))+((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)/d;
			C(3, 5) = (rCG(1, 2)+rCG(1, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))+((rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)/d;
			C(4, 0) = ((rCG(1, 2)+rCG(1, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*-2.0)/d;
			C(4, 1) = ((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(4, 2) = ((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(4, 3) = (rCG(0, 2)+rCG(0, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))+((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)/d;
			C(4, 4) = -(rCG(0, 0)-(pow(rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0,2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0+rCG(0, 0)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/d;
			C(4, 5) = (rCG(0, 1)+rCG(0, 1)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))+((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)/d;
			C(5, 0) = ((rCG(1, 1)*rCG(2, 2)-rCG(1, 2)*rCG(1, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(5, 1) = ((rCG(0, 2)+rCG(0, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))-((rCG(0, 0)*rCG(2, 2)-rCG(0, 2)*rCG(0, 2))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)*-2.0)/d;
			C(5, 2) = ((rCG(0, 0)*rCG(1, 1)-rCG(0, 1)*rCG(0, 1))*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/d;
			C(5, 3) = (rCG(1, 2)+rCG(1, 2)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))+((rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*(rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 1)*rCG(2, 2)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)/d;
			C(5, 4) = (rCG(0, 1)+rCG(0, 1)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2))+((rCG(0, 1)*rCG(0, 2)*2.0-rCG(0, 0)*rCG(1, 2)*2.0)*(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0)/d;
			C(5, 5) = -(rCG(1, 1)-(pow(rCG(0, 1)*rCG(1, 2)*2.0-rCG(0, 2)*rCG(1, 1)*2.0,2.0)*1.0/pow(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2),2.0))/2.0+rCG(1, 1)/(rCG(0, 0)*(rCG(1, 2)*rCG(1, 2))+(rCG(0, 2)*rCG(0, 2))*rCG(1, 1)+(rCG(0, 1)*rCG(0, 1))*rCG(2, 2)-rCG(0, 1)*rCG(0, 2)*rCG(1, 2)*2.0-rCG(0, 0)*rCG(1, 1)*rCG(2, 2)))/d;
			break;
		case HyperElasticPropertiesConfiguration::MODEL::BLATZ_KO_FOAM:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::GENT:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::OGDEN_1:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::OGDEN_2:
			break;
		case HyperElasticPropertiesConfiguration::MODEL::OGDEN_3:
			break;
		}
	}
}

// source dX, dY, dZ

// target::
// dX   0   0
//  0  dY   0
//  0   0  dZ
// dY  dX   0
//  0  dZ  dY
// dZ   0  dX
static void distribute6x3(double *target, const double *source, size_t rows, size_t columns)
{
	const double *dNDx = source;
	const double *dNDy = source + columns;
	const double *dNDz = source + 2 * columns;

	memcpy(target                                   , dNDx, sizeof(double) * columns);
	memcpy(target + 3 * rows * columns +     columns, dNDx, sizeof(double) * columns);
	memcpy(target + 5 * rows * columns + 2 * columns, dNDx, sizeof(double) * columns);

	memcpy(target + 1 * rows * columns +     columns, dNDy, sizeof(double) * columns);
	memcpy(target + 3 * rows * columns              , dNDy, sizeof(double) * columns);
	memcpy(target + 4 * rows * columns + 2 * columns, dNDy, sizeof(double) * columns);

	memcpy(target + 2 * rows * columns + 2 * columns, dNDz, sizeof(double) * columns);
	memcpy(target + 4 * rows * columns +     columns, dNDz, sizeof(double) * columns);
	memcpy(target + 5 * rows * columns              , dNDz, sizeof(double) * columns);
}

// source dX, dY, dZ

// target::
// dX   0   0
// dY   0   0
// dZ   0   0
//  0  dX   0
//  0  dY   0
//  0  dZ   0
//  0   0  dX
//  0   0  dY
//  0   0  dZ
static void distribute9x3(double *target, const double *source, size_t rows, size_t columns)
{
	const double *dNDx = source;
	const double *dNDy = source + columns;
	const double *dNDz = source + 2 * columns;

	memcpy(target                                   , dNDx, sizeof(double) * columns);
	memcpy(target + 3 * rows * columns + 1 * columns, dNDx, sizeof(double) * columns);
	memcpy(target + 6 * rows * columns + 2 * columns, dNDx, sizeof(double) * columns);

	memcpy(target + 1 * rows * columns              , dNDy, sizeof(double) * columns);
	memcpy(target + 4 * rows * columns + 1 * columns, dNDy, sizeof(double) * columns);
	memcpy(target + 7 * rows * columns + 2 * columns, dNDy, sizeof(double) * columns);

	memcpy(target + 2 * rows * columns              , dNDz, sizeof(double) * columns);
	memcpy(target + 5 * rows * columns + 1 * columns, dNDz, sizeof(double) * columns);
	memcpy(target + 8 * rows * columns + 2 * columns, dNDz, sizeof(double) * columns);
}

void StructuralMechanics3DKernel::processElement(const Builder &builder, const ElasticityElementIterator &iterator, InstanceFiller &filler) const
{
	iterator.corotating = NULL;
	iterator.fixed = NULL;
	int rsize = info::mesh->elements->regions->edataSize();
	for (size_t r = 0; r < info::mesh->elementsRegions.size(); r++) {
		esint maskOffset = r / (8 * sizeof(esint));
		esint bit = (esint)1 << (r % (8 * sizeof(esint)));
		if (info::mesh->elements->regions->datatarray()[iterator.offset * rsize + maskOffset] & bit) {
			switch (iterator.configuration.rotor_dynamics.type) {
			case RotorDynamicsConfiguration::TYPE::FIXED: {
				auto def = iterator.configuration.rotor_dynamics.fixed.rotors_definitions.find(info::mesh->elementsRegions[r]->name);
				if (def != iterator.configuration.rotor_dynamics.fixed.rotors_definitions.end()) {
					iterator.fixed = &def->second;
				}
			} break;
			case RotorDynamicsConfiguration::TYPE::COROTATING: {
				auto def = iterator.configuration.rotor_dynamics.corotating.rotors_definitions.find(info::mesh->elementsRegions[r]->name);
				if (def != iterator.configuration.rotor_dynamics.corotating.rotors_definitions.end()) {
					iterator.corotating = &def->second;
					iterator.rotationAxis = &iterator.configuration.rotor_dynamics.corotating.rotation_axis;
				}
			} break;
			}
		}
	}

	bool harmonic = step::step.type == step::TYPE::FREQUENCY || step::step.type == step::TYPE::FTT;

	int size = iterator.element->nodes;
	filler.DOFs = 3 * size;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &NNN = *(iterator.element->NNN);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	MatrixDense C(6, 6), XYZ(1, 3), initCoordinates(size, 3), coordinates(size, 3), F(3, 3), J, JC, invJ(3, 3), dND, B, CEp, CB, precision, rhsT, SGL, CBL, disp(3 * size, 1), prestress(3 * size, 3 * size);
	MatrixDense eHat(3, 3), eVec(6, 1), sVec(6, 1), S(9, 9), BL(6, 3 * size), GL(9, 3 * size);
	MatrixDense K(size, 36), TE(size, 3), inertia(size, 3), dens(size, 1);
	MatrixDense gpK(size, 36), gpTE(1, 3), gpInertia(1, 3), gpDens(1, 1);
	MatrixDense rotation(3, 3), spin(3, 3), Ks, omegaN;
	MatrixDense uc(3 * size, 1), us(3 * size, 1), uB(6, 1);
	Point fixedOmega, fixedP;
	MatrixDense tx(1, 3), ty(1, 3), tz(1, 3), G(3 * size, 3 * size), Nxr(3, 1), x, fixedR(3, 1), tztyNxr(3, 1), tytzNxr(3, 1), BB(1, size), Bx(3 * size, 3), Bxt1(3 * size, 1), Bxt2(3 * size, 1), Bxttt(3 * size, 3);
	double detJ, te;

	Point orientation;
	if (iterator.material->linear_elastic_properties.orientation) {
		orientation.x = this->orientation->data[3 * iterator.offset + 0];
		orientation.y = this->orientation->data[3 * iterator.offset + 1];
		orientation.z = this->orientation->data[3 * iterator.offset + 2];
	}

	for (int n = 0; n < size; n++) {
		Evaluator::Params params;
		params.coords(3, iterator.coordinates.data + 3 * n);
		params.temp(iterator.temperature.data + n);
		inertia(n, 0) = iterator.acceleration.data[3 * n + 0];
		inertia(n, 1) = iterator.acceleration.data[3 * n + 1];
		inertia(n, 2) = iterator.acceleration.data[3 * n + 2];
		initCoordinates(n, 0) = iterator.coordinates.data[3 * n + 0];
		initCoordinates(n, 1) = iterator.coordinates.data[3 * n + 1];
		initCoordinates(n, 2) = iterator.coordinates.data[3 * n + 2];
		coordinates(n, 0) = initCoordinates(n, 0) + iterator.displacement.data[3 * n + 0];
		coordinates(n, 1) = initCoordinates(n, 1) + iterator.displacement.data[3 * n + 1];
		coordinates(n, 2) = initCoordinates(n, 2) + iterator.displacement.data[3 * n + 2];
		disp(n + 0 * size, 0) = iterator.displacement.data[3 * n + 0];
		disp(n + 1 * size, 0) = iterator.displacement.data[3 * n + 1];
		disp(n + 2 * size, 0) = iterator.displacement.data[3 * n + 2];
		dens(n, 0) = iterator.material->density.evaluator->eval(params);

		if (harmonic) {
			uc(0 * size + n, 0) = iterator.cos.data[3 * n + 0];
			uc(1 * size + n, 0) = iterator.cos.data[3 * n + 1];
			uc(2 * size + n, 0) = iterator.cos.data[3 * n + 2];

			us(0 * size + n, 0) = iterator.sin.data[3 * n + 0];
			us(1 * size + n, 0) = iterator.sin.data[3 * n + 1];
			us(2 * size + n, 0) = iterator.sin.data[3 * n + 2];
		}

		switch (iterator.material->thermal_expansion.model) {
		case ThermalExpansionConfiguration::MODEL::ISOTROPIC:
			te = iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->eval(params);
			TE(n, 0) = TE(n, 1) = TE(n, 2) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			break;
		case ThermalExpansionConfiguration::MODEL::ORTHOTROPIC:
			te = iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->eval(params);
			TE(n, 0) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			te = iterator.material->thermal_expansion.thermal_expansion.get(1, 1).evaluator->eval(params);
			TE(n, 1) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			te = iterator.material->thermal_expansion.thermal_expansion.get(2, 2).evaluator->eval(params);
			TE(n, 2) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			break;
		}

		switch (iterator.material->material_model) {
		case MaterialBaseConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC:
			assembleLinearElasticMaterialMatrix(n, iterator.coordinates.data, iterator.material, step::time.current, iterator.temperature.data[n], K, orientation);
			break;
		case MaterialBaseConfiguration::MATERIAL_MODEL::HYPER_ELASTIC:
//			assembleHyperElasticMaterialMatrix(n, iterator.coordinates, iterator.material, step::time.current, iterator.temperature[n], K);
			break;
		}
	}

	x.nrows = 3 * size;
	x.ncols = 1;
	x.vals = initCoordinates.vals;

	if (iterator.fixed) {
		fixedOmega.x = iterator.fixed->rotation_axis.orientation.x.evaluator->eval(Evaluator::Params());
		fixedOmega.y = iterator.fixed->rotation_axis.orientation.y.evaluator->eval(Evaluator::Params());
		fixedOmega.z = iterator.fixed->rotation_axis.orientation.z.evaluator->eval(Evaluator::Params());
		fixedOmega.normalize();
		fixedP.x = fixedOmega.x;
		fixedP.y = fixedOmega.y;
		fixedR[0][0] = iterator.fixed->rotation_axis.center.x.evaluator->eval(Evaluator::Params());
		fixedR[1][0] = iterator.fixed->rotation_axis.center.y.evaluator->eval(Evaluator::Params());
		fixedR[2][0] = iterator.fixed->rotation_axis.center.z.evaluator->eval(Evaluator::Params());
		double cosphi = 1, sinphi = 0, costheta = 0, sintheta = fixedOmega.z < 0 ? -1 : 1;
		if (fixedP.norm() > 1e-9) {
			cosphi = fixedP.x / fixedP.norm();
			sinphi = fixedOmega.y < 0 ? - std::sqrt(1 - cosphi * cosphi) : std::sqrt(1 - cosphi * cosphi);
			costheta = Point(cosphi, sinphi, 0) * fixedOmega;
			sintheta = fixedOmega.z < 0 ? - std::sqrt(1 - costheta * costheta) : std::sqrt(1 - costheta * costheta);
		}
		tx[0][0] = cosphi * costheta;
		tx[0][1] = sinphi * costheta;
		tx[0][2] = sintheta;

		ty[0][0] = -sinphi;
		ty[0][1] = cosphi;
		ty[0][2] = 0;

		tz[0][0] = -cosphi * sintheta;
		tz[0][1] = -sinphi * sintheta;
		tz[0][2] = costheta;
	}

	if (iterator.corotating && iterator.corotating->coriolis_effect) {
		double ox = step::frequency.angular / iterator.corotating->rotation.frequency_ratio * iterator.rotationAxis->orientation.x.evaluator->eval(Evaluator::Params());
		double oy = step::frequency.angular / iterator.corotating->rotation.frequency_ratio * iterator.rotationAxis->orientation.y.evaluator->eval(Evaluator::Params());
		double oz = step::frequency.angular / iterator.corotating->rotation.frequency_ratio * iterator.rotationAxis->orientation.z.evaluator->eval(Evaluator::Params());
		rotation.fill(0);
		rotation(0, 1) = -oz;
		rotation(0, 2) =  oy;
		rotation(1, 0) =  oz;
		rotation(1, 2) = -ox;
		rotation(2, 0) = -oy;
		rotation(2, 1) =  ox;
//		std::cout << rotation;
	}

	if (iterator.corotating && iterator.corotating->spin_softening) {
		double ox = step::frequency.angular / iterator.corotating->rotation.frequency_ratio * iterator.rotationAxis->orientation.x.evaluator->eval(Evaluator::Params());
		double oy = step::frequency.angular / iterator.corotating->rotation.frequency_ratio * iterator.rotationAxis->orientation.y.evaluator->eval(Evaluator::Params());
		double oz = step::frequency.angular / iterator.corotating->rotation.frequency_ratio * iterator.rotationAxis->orientation.z.evaluator->eval(Evaluator::Params());
		spin(0, 0) = -(oy * oy + oz * oz);
		spin(1, 1) = -(ox * ox + oz * oz);
		spin(2, 2) = -(ox * ox + oy * oy);
		spin(0, 1) = ox * oy;
		spin(0, 2) = ox * oz;
		spin(1, 0) = ox * oy;
		spin(1, 2) = oy * oz;
		spin(2, 0) = ox * oz;
		spin(2, 1) = oy * oz;
	}

	if (builder.matrices & (Builder::Request::K | Builder::Request::R)) {
		filler.Ke.resize(3 * size, 3 * size);
		filler.Ke.fill(0);
	}
	if (builder.matrices & Builder::Request::C) {
		filler.Ce.resize(3 * size, 3 * size);
		filler.Ce.fill(0);
		filler.CMe.resize(3 * size, 3 * size);
		filler.CMe.fill(0);
	}
	if (builder.matrices & (Builder::Request::M | Builder::Request::R)) {
		filler.Me.resize(3 * size, 3 * size);
		filler.Me.fill(0);
	}
	if (builder.matrices & Builder::Request::R) {
		filler.Re.resize(3 * size);
		filler.Re.fill(0);
	}
	if (builder.matrices & Builder::Request::f) {
		filler.Fe.resize(3 * size);
		filler.Fe.fill(0);
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], initCoordinates);

		detJ = MATH::determinant3x3(J.vals);
		if (detJ < 0) { ++filler.invalid; detJ = -detJ; }
		MATH::Dense3x3inverse(J.vals, invJ.vals, detJ);

		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);
		gpDens.multiply(N[gp], dens);

		if (builder.matrices & Builder::Request::f) {
			gpTE.multiply(N[gp], TE);
			gpInertia.multiply(N[gp], inertia);
			XYZ.multiply(N[gp], coordinates);
		}

		if (builder.matrices & Builder::Request::C) {
			if (iterator.corotating && iterator.corotating->coriolis_effect) {
				omegaN.multiply(rotation, NNN[gp]);
				filler.Ce.multiply(NNN[gp], omegaN, 2 * gpDens(0, 0) * detJ * weighFactor[gp], 1, true);
			}
		}

		if (iterator.corotating && iterator.corotating->spin_softening) {
			omegaN.multiply(spin, NNN[gp]);
			Ks.multiply(NNN[gp], omegaN, gpDens(0, 0) * detJ * weighFactor[gp], 1, true);
		}

		if (step::step.type != step::TYPE::FTT && (builder.matrices & (Builder::Request::M | Builder::Request::R))) {
			filler.Me.multiply(NNN[gp], NNN[gp], gpDens(0, 0) * detJ * weighFactor[gp], 1, true);
		}

		if (iterator.material->material_model == MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC) {
			size_t k = 0;
			for (size_t i = 0; i < 6; i++) {
				C(i, i) = gpK(0, k++);
			}
			for (size_t i = 0; i < 6; i++) {
				for (size_t j = i + 1; j < 6; j++) {
					C(i, j) = gpK(0, k++);
				}
			}
			for (size_t i = 0; i < 6; i++) {
				for (size_t j = 0; j < i; j++) {
					C(i, j) = gpK(0, k++);
				}
			}
		}

		if (builder.prestress) {
			B.resize(C.nrows, 3 * size);
			B.fill(0);
			distribute6x3(B.vals, dND.vals, dND.nrows, dND.ncols);
			// 6x24 * 24x1
			uB.multiply(B, disp);
			uB(0, 0) -= gpTE(0, 0);
			uB(1, 0) -= gpTE(0, 1);
			uB(2, 0) -= gpTE(0, 2);

			MatrixDense stress(3, 3), Sigma, stress88(size, size), stressXX;
			Sigma.multiply(C, uB);
			stress[0][0] = Sigma[0][0];
			stress[1][1] = Sigma[1][0];
			stress[2][2] = Sigma[2][0];
			stress[0][1] = Sigma[3][0];
			stress[1][0] = Sigma[3][0];
			stress[1][2] = Sigma[4][0];
			stress[2][1] = Sigma[4][0];
			stress[0][2] = Sigma[5][0];
			stress[2][0] = Sigma[5][0];

			stressXX.multiply(stress, dND);
			stress88.multiply(dND, stressXX, 1, 0, true, false);
			for (int i = 0; i < size; ++i) {
				for (int j = 0; j < size; ++j) {
					prestress[i][j] = stress88[i][j];
					prestress[i + size][j + size] = stress88[i][j];
					prestress[i + 2 * size][j + 2 * size] = stress88[i][j];
				}
			}
		}

		if (iterator.material->material_model == MaterialConfiguration::MATERIAL_MODEL::HYPER_ELASTIC && step::isInitial()) {
			assembleHyperElasticMaterialMatrix(iterator.material, F, C, S);
		}

		if (iterator.largeDisplacement && !step::isInitial()) {
			JC.multiply(dN[gp], coordinates);
			F.multiply(JC, invJ, 1, 0, true);

			if (iterator.material->material_model == MaterialConfiguration::MATERIAL_MODEL::HYPER_ELASTIC) {
				assembleHyperElasticMaterialMatrix(iterator.material, F, C, S);
				for (int i = 1; i < 3; i++) {
					S(0 + i * 3, 0 + i * 3) = S(0, 0);
					S(1 + i * 3, 1 + i * 3) = S(1, 1);
					S(2 + i * 3, 2 + i * 3) = S(2, 2);
					S(0 + i * 3, 1 + i * 3) = S(0, 1);
					S(1 + i * 3, 0 + i * 3) = S(1, 0);
					S(0 + i * 3, 2 + i * 3) = S(0, 2);
					S(2 + i * 3, 0 + i * 3) = S(2, 0);
					S(1 + i * 3, 2 + i * 3) = S(1, 2);
					S(2 + i * 3, 1 + i * 3) = S(2, 1);
				}
				sVec(0, 0) = S(0, 0);
				sVec(1, 0) = S(1, 1);
				sVec(2, 0) = S(2, 2);
				sVec(3, 0) = S(0, 1);
				sVec(4, 0) = S(1, 2);
				sVec(5, 0) = S(0, 2);
			}

			if (iterator.material->material_model == MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC) {
				eHat.multiply(F, F, .5, 0, true, false);
				eHat(0, 0) -= .5;
				eHat(1, 1) -= .5;
				eHat(2, 2) -= .5;
				eVec(0, 0) = eHat(0, 0);
				eVec(1, 0) = eHat(1, 1);
				eVec(2, 0) = eHat(2, 2);
				eVec(3, 0) = 2 * eHat(0, 1);
				eVec(4, 0) = 2 * eHat(1, 2);
				eVec(5, 0) = 2 * eHat(0, 2);
				sVec.multiply(C, eVec);
				for (int i = 0; i < 3; i++) {
					S(0 + i * 3, 0 + i * 3) = sVec(0, 0);
					S(1 + i * 3, 1 + i * 3) = sVec(1, 0);
					S(2 + i * 3, 2 + i * 3) = sVec(2, 0);
					S(0 + i * 3, 1 + i * 3) = sVec(3, 0);
					S(1 + i * 3, 0 + i * 3) = sVec(3, 0);
					S(0 + i * 3, 2 + i * 3) = sVec(5, 0);
					S(2 + i * 3, 0 + i * 3) = sVec(5, 0);
					S(1 + i * 3, 2 + i * 3) = sVec(4, 0);
					S(2 + i * 3, 1 + i * 3) = sVec(4, 0);
				}
			}

			distribute9x3(GL.vals, dND.vals, dND.nrows, dND.ncols);
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < 3; j++) {
					BL(0, i + j * size) = F(j, 0) * dND(0, i);
					BL(1, i + j * size) = F(j, 1) * dND(1, i);
					BL(2, i + j * size) = F(j, 2) * dND(2, i);
					BL(3, i + j * size) = F(j, 0) * dND(1, i) + F(j, 1) * dND(0, i);
					BL(4, i + j * size) = F(j, 1) * dND(2, i) + F(j, 2) * dND(1, i);
					BL(5, i + j * size) = F(j, 0) * dND(2, i) + F(j, 2) * dND(0, i);
				}
			}

			if (builder.matrices & Builder::Request::K) {
				CBL.multiply(C, BL);
				filler.Ke.multiply(BL, CBL, detJ * weighFactor[gp], 1, true);
				SGL.multiply(S, GL);
				filler.Ke.multiply(GL, SGL, detJ * weighFactor[gp], 1, true);
			}
			if (builder.matrices & Builder::Request::R) {
				filler.Re.multiply(BL, sVec, detJ * weighFactor[gp], 1, true);
			}
		} else {
			B.resize(C.nrows, 3 * size);
			B.fill(0);
			distribute6x3(B.vals, dND.vals, dND.nrows, dND.ncols);
			if (step::step.type != step::TYPE::FTT && (builder.matrices & (Builder::Request::K | Builder::Request::R))) {
				CB.multiply(C, B);
				filler.Ke.multiply(B, CB, detJ * weighFactor[gp], 1, true);
				if (iterator.fixed) {
					B.multiply(tx, dND);
					for (esint c = 0; c < B.ncols; c++) {
						Bx[0 * B.ncols + c][0] = B[0][c];
						Bx[1 * B.ncols + c][1] = B[0][c];
						Bx[2 * B.ncols + c][2] = B[0][c];
					}
					Nxr.multiply(NNN[gp], x);
					Nxr.add(-1, &fixedR);
					double tyNxr = ty[0][0] * Nxr[0][0] + ty[0][1] * Nxr[1][0] + ty[0][2] * Nxr[2][0];
					double tzNxr = tz[0][0] * Nxr[0][0] + tz[0][1] * Nxr[1][0] + tz[0][2] * Nxr[2][0];
					Bxt1.multiply(Bx, tz, tyNxr, 0, false, true);
					Bxt2.multiply(Bx, ty, tzNxr, 0, false, true);
					Bxt1.add(-1, &Bxt2);
					Bxttt.multiply(Bxt1, tx);
					G.multiply(Bxttt, NNN[gp], weighFactor[gp] * detJ * gpDens[0][0], 1);
				}
			}

			if (step::step.type != step::TYPE::FTT && (builder.matrices & Builder::Request::f)) {
				precision.resize(C.nrows, 1);
				precision(0, 0) = gpTE(0, 0);
				precision(1, 0) = gpTE(0, 1);
				precision(2, 0) = gpTE(0, 2);
				precision(3, 0) = precision(4, 0) = precision(5, 0) = 0;

				CEp.multiply(C, precision);
				rhsT.multiply(B, CEp, detJ * weighFactor[gp], 0, true, false);
				for (esint i = 0; i < 3 * size; i++) {
					filler.Fe[0][i] += rhsT(i, 0);
				}
			}
		}

		for (esint i = 0; i < 3 * size; i++) {
			filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % size) * gpInertia(0, i / size);
			if (iterator.angularVelocity.data[0]) {
				filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % size) * (i / size == 0 ? 0 : XYZ(0, i / size)) * pow(iterator.angularVelocity.data[0], 2);
			}
			if (iterator.angularVelocity.data[1]) {
				filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % size) * (i / size == 1 ? 0 : XYZ(0, i / size)) * pow(iterator.angularVelocity.data[1], 2);
			}
			if (iterator.angularVelocity.data[2]) {
				filler.Fe[0][i] += gpDens(0, 0) * detJ * weighFactor[gp] * N[gp](0, i % size) * (i / size == 2 ? 0 : XYZ(0, i / size)) * pow(iterator.angularVelocity.data[2], 2);
			}
		}
	}

	if (builder.matrices & Builder::Request::C) {
		if (iterator.fixed) {
			filler.Ce.add(step::frequency.angular * iterator.fixed->rotation.frequency_ratio, &G);
			G.transpose();
			filler.Ce.add(-step::frequency.angular * iterator.fixed->rotation.frequency_ratio, &G);
		}
		if (builder.rayleighDamping) {
			double stiffDamping = builder.stiffnessDamping + builder.structuralDampingCoefficient / step::frequency.angular;
			filler.Ce.add(step::frequency.angular * stiffDamping, &filler.Ke);
			filler.Ce.add(step::frequency.angular * builder.massDamping, &filler.Me);
		}
	}

	if (builder.prestress && (builder.matrices & Builder::Request::K)) {
		filler.Ke.add(1, &prestress);
	}

	if (iterator.corotating && iterator.corotating->spin_softening) {
		filler.Ke.add(-1, &Ks);
	}

	if (iterator.massStabilization) {
		auto dual = info::mesh->elements->faceNeighbors->begin() + iterator.offset;
		auto fpointer = iterator.element->facepointers->begin();
		auto fnodes = iterator.element->faces->begin();
		auto doffset = std::lower_bound(info::mesh->domains->elements.begin(), info::mesh->domains->elements.end(), iterator.offset + 1) - info::mesh->domains->elements.begin() - 1;
		auto lower = info::mesh->domains->elements[doffset] + info::mesh->elements->distribution.process.offset;
		auto upper = info::mesh->domains->elements[doffset + 1] + info::mesh->elements->distribution.process.offset;
		for (auto neigh = dual->begin(); neigh != dual->end(); ++neigh, ++fpointer, ++fnodes) {
			if (*neigh != -1) {
				if (*neigh < lower || upper <= *neigh) {
					double sign = 1;
					if (*neigh < lower) {
						sign = -1;
					}
					MatrixDense cc(fnodes->size(), 3);
					for (size_t n = 0; n < fnodes->size(); ++n) {
						cc(n, 0) = iterator.coordinates.data[3 * fnodes->at(n) + 0];
						cc(n, 1) = iterator.coordinates.data[3 * fnodes->at(n) + 1];
						cc(n, 2) = iterator.coordinates.data[3 * fnodes->at(n) + 2];
					}
					MatrixDense fdND;
					for (size_t fgp = 0; fgp < fpointer->at(0)->dN->size(); ++fgp) {
						fdND.multiply(fpointer->at(0)->dN->at(fgp), coordinates);
						Point v2(fdND(0, 0), fdND(0, 1), fdND(0, 2));
						Point v1(fdND(1, 0), fdND(1, 1), fdND(1, 2));
						Point va = Point::cross(v1, v2);
						double fJ = va.norm();
						filler.CMe.multiply(fpointer->at(0)->NNN->at(fgp), fpointer->at(0)->NNN->at(fgp), sign * fpointer->at(0)->weighFactor->at(fgp) * fJ, 1, true);
					}
				}
			}
		}
	}

	if (!iterator.largeDisplacement && step::step.type != step::TYPE::FTT && (builder.matrices & Builder::Request::R)) {
		if (harmonic) {
			filler.Re[0].multiply(filler.Ke, uc, 1.0);
			filler.Re[1].multiply(filler.Ke, us, 1.0);
			filler.Re[0].multiply(filler.Me, uc, -step::frequency.angular * step::frequency.angular, 1.0);
			filler.Re[1].multiply(filler.Me, us, -step::frequency.angular * step::frequency.angular, 1.0);

			if (builder.matrices & Builder::Request::C) {
				filler.Re[0].multiply(filler.Ce, us, -1, 1.0);
				filler.Re[1].multiply(filler.Ce, uc,  1, 1.0);
			} else {
				filler.Re[0].multiply(filler.Ke, us, -step::frequency.angular * (builder.stiffnessDamping + builder.structuralDampingCoefficient / step::frequency.angular), 1.0);
				filler.Re[0].multiply(filler.Me, us, -step::frequency.angular * builder.massDamping, 1.0);

				filler.Re[1].multiply(filler.Ke, uc,  step::frequency.angular * (builder.stiffnessDamping + builder.structuralDampingCoefficient / step::frequency.angular), 1.0);
				filler.Re[1].multiply(filler.Me, uc,  step::frequency.angular * builder.massDamping, 1.0);
			}

			if (!(builder.matrices & Builder::Request::M)) {
				filler.Me.resize(0, 0);
			}
			if (!(builder.matrices & Builder::Request::K)) {
				filler.Ke.resize(0, 0);
			}
		}
	}

	filler.insertK = builder.matrices & Builder::Request::K;
	filler.insertM = builder.matrices & Builder::Request::M;
	filler.insertC = builder.matrices & Builder::Request::C;
	filler.insertR = builder.matrices & Builder::Request::R;
	filler.insertF = builder.matrices & Builder::Request::f;
}

void StructuralMechanics3DKernel::processFace(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	filler.insertK = filler.insertF = false;
	filler.DOFs = 3 * iterator.element->nodes;
	if (iterator.normalPressure.data == NULL) {
		return;
	}
	if (!(builder.matrices & (Builder::Request::K | Builder::Request::f))) {
		return;
	}

	esint size = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);
	const std::vector<double> &weighFactor = *(iterator.element->weighFactor);

	MatrixDense coordinates(size, 3), dND(1, 3), P(size, 1), normal(1, 3);
	MatrixDense gpP(1, 1), gpQ(1, 3);

	filler.DOFs = 3 * size;
	if ((filler.insertF = (builder.matrices & Builder::Request::f))) {
		filler.Fe.resize(filler.DOFs);
		filler.Fe.fill(0);
	}

	for (esint n = 0; n < size; n++) {
		coordinates(n, 0) = iterator.coordinates.data[3 * n + 0];
		coordinates(n, 1) = iterator.coordinates.data[3 * n + 1];
		coordinates(n, 2) = iterator.coordinates.data[3 * n + 2];
		P(n, 0) += iterator.normalPressure.data[n];
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		dND.multiply(dN[gp], coordinates);
		Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
		Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
		Point va = Point::cross(v1, v2);
		// e->rotateOutside(e->parentElements()[0], _mesh->coordinates(), va);
		double J = va.norm();
		normal(0, 0) = va.x / va.norm();
		normal(0, 1) = va.y / va.norm();
		normal(0, 2) = va.z / va.norm();

		gpP.multiply(N[gp], P);
		gpQ.multiply(normal, gpP, 1, 0, true);

		for (esint i = 0; i < filler.DOFs; i++) {
			filler.Fe[0][i] += J * weighFactor[gp] * N[gp](0, i % size) * gpQ(0, i / size);
		}
	}
}

void StructuralMechanics3DKernel::processEdge(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	filler.insertK = filler.insertF = false;
	filler.DOFs = 3 * iterator.element->nodes;
}

void StructuralMechanics3DKernel::processNode(const Builder &builder, const ElasticityBoundaryIterator &iterator, InstanceFiller &filler) const
{
	filler.insertK = filler.insertF = false;
	filler.DOFs = 3;
	if (iterator.force.data == NULL && iterator.harmonicForce.magnitude.data == NULL && iterator.rotatingForce.mass.data == NULL) {
		return;
	}
	if (!(builder.matrices & Builder::Request::f)) {
		return;
	}
	if (step::step.type == step::TYPE::FTT) {
		return;
	}

	if ((filler.insertF = (builder.matrices & Builder::Request::f))) {
		filler.Fe.resize(3);
		filler.Fe.fill(0);
	}

	if (iterator.force.data) {
		filler.Fe[0][0] = iterator.force.data[0];
		filler.Fe[0][1] = iterator.force.data[1];
		filler.Fe[0][2] = iterator.force.data[2];
	}

	if (iterator.harmonicForce.magnitude.data != NULL) {
		filler.Fe[0][0] = iterator.harmonicForce.magnitude.data[0] * std::cos(iterator.harmonicForce.phase.data[0] * M_PI / 180);
		filler.Fe[0][1] = iterator.harmonicForce.magnitude.data[1] * std::cos(iterator.harmonicForce.phase.data[1] * M_PI / 180);
		filler.Fe[0][2] = iterator.harmonicForce.magnitude.data[2] * std::cos(iterator.harmonicForce.phase.data[2] * M_PI / 180);
		filler.Fe[1][0] = iterator.harmonicForce.magnitude.data[0] * std::sin(iterator.harmonicForce.phase.data[0] * M_PI / 180);
		filler.Fe[1][1] = iterator.harmonicForce.magnitude.data[1] * std::sin(iterator.harmonicForce.phase.data[1] * M_PI / 180);
		filler.Fe[1][2] = iterator.harmonicForce.magnitude.data[2] * std::sin(iterator.harmonicForce.phase.data[2] * M_PI / 180);
	}

	if (iterator.rotatingForce.mass.data != NULL) {
		double value = step::frequency.angular * step::frequency.angular * iterator.rotatingForce.mass.data[0] * iterator.rotatingForce.radius.data[0];
		filler.Fe[0][0] =   value * std::cos(iterator.rotatingForce.phaseAngle.data[0] * M_PI / 180);
		filler.Fe[0][1] = - value * std::sin(iterator.rotatingForce.phaseAngle.data[0] * M_PI / 180);
		filler.Fe[1][0] = - value * std::sin(iterator.rotatingForce.phaseAngle.data[0] * M_PI / 180);
		filler.Fe[1][1] = - value * std::cos(iterator.rotatingForce.phaseAngle.data[0] * M_PI / 180);
	}
}

void StructuralMechanics3DKernel::elementSolution(ElasticityElementIterator &iterator)
{
	if (!info::ecf->output.results_selection.stress || step::step.type == step::TYPE::FREQUENCY) {
		return;
	}

	Point orientation;
	if (iterator.material->linear_elastic_properties.orientation) {
		orientation.x = this->orientation->data[3 * iterator.offset + 0];
		orientation.x = this->orientation->data[3 * iterator.offset + 1];
		orientation.x = this->orientation->data[3 * iterator.offset + 2];
	}

	int size = iterator.element->nodes;

	const std::vector<MatrixDense> &N = *(iterator.element->N);
	const std::vector<MatrixDense> &dN = *(iterator.element->dN);

	MatrixDense C(6, 6), initCoordinates(size, 3), disp(3 * size, 1), J, invJ(3, 3), dND, B, F(3, 3), S(3, 3), SF(3, 3), Sigma(6, 1);
	MatrixDense K(size, 36), TE(size, 3), uB(6, 1);
	MatrixDense gpK(size, 36), gpTE(1, 3), gpInertia(1, 3), gpDens(1, 1);
	double detJ, detF, te;

	for (int n = 0; n < size; n++) {
		Evaluator::Params params;
		params.coords(3, iterator.coordinates.data + 3 * n);
		params.temp(iterator.temperature.data + n);
		initCoordinates(n, 0) = iterator.coordinates.data[3 * n + 0];
		initCoordinates(n, 1) = iterator.coordinates.data[3 * n + 1];
		initCoordinates(n, 2) = iterator.coordinates.data[3 * n + 2];
		disp(n + 0 * size, 0) = iterator.displacement.data[3 * n + 0];
		disp(n + 1 * size, 0) = iterator.displacement.data[3 * n + 1];
		disp(n + 2 * size, 0) = iterator.displacement.data[3 * n + 2];

		switch (iterator.material->thermal_expansion.model) {
		case ThermalExpansionConfiguration::MODEL::ISOTROPIC:
			te = iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->eval(params);
			TE(n, 0) = TE(n, 1) = TE(n, 2) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			break;
		case ThermalExpansionConfiguration::MODEL::ORTHOTROPIC:
			te = iterator.material->thermal_expansion.thermal_expansion.get(0, 0).evaluator->eval(params);
			TE(n, 0) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			te = iterator.material->thermal_expansion.thermal_expansion.get(1, 1).evaluator->eval(params);
			TE(n, 1) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			te = iterator.material->thermal_expansion.thermal_expansion.get(2, 2).evaluator->eval(params);
			TE(n, 2) = (iterator.temperature.data[n] - iterator.initialTemperature.data[n]) * te;
			break;
		}

		switch (iterator.material->material_model) {
		case MaterialBaseConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC:
			assembleLinearElasticMaterialMatrix(n, iterator.coordinates.data, iterator.material, step::time.current, iterator.temperature.data[n], K, orientation);
			break;
		case MaterialBaseConfiguration::MATERIAL_MODEL::HYPER_ELASTIC:
//			assembleHyperElasticMaterialMatrix(n, iterator.coordinates, iterator.material, step::time.current, iterator.temperature[n], K);
			break;
		}
	}

	for (size_t gp = 0; gp < N.size(); gp++) {
		J.multiply(dN[gp], initCoordinates);
		detJ = MATH::determinant3x3(J.vals);
		if (detJ < 0) { detJ = -detJ; }
		MATH::Dense3x3inverse(J.vals, invJ.vals, detJ);

		gpK.multiply(N[gp], K);
		dND.multiply(invJ, dN[gp]);
		gpTE.multiply(N[gp], TE);

		if (iterator.material->material_model == MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC) {
			size_t k = 0;
			for (size_t i = 0; i < 6; i++) {
				C(i, i) = gpK(0, k++);
			}
			for (size_t i = 0; i < 6; i++) {
				for (size_t j = i + 1; j < 6; j++) {
					C(i, j) = gpK(0, k++);
				}
			}
			for (size_t i = 0; i < 6; i++) {
				for (size_t j = 0; j < i; j++) {
					C(i, j) = gpK(0, k++);
				}
			}
		}

		if (iterator.material->material_model == MaterialConfiguration::MATERIAL_MODEL::HYPER_ELASTIC) {
			assembleHyperElasticMaterialMatrix(iterator.material, F, C, S);
			detF = MATH::determinant3x3(F.vals);
			SF.multiply(S, F, 1, 0, false, true);

			Sigma.multiply(F, SF, 1 / detF, 1);
		}
		if (iterator.material->material_model == MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC) {
			B.resize(C.nrows, 3 * size);
			B.fill(0);
			distribute6x3(B.vals, dND.vals, dND.nrows, dND.ncols);
			// 6x24 * 24x1
			uB.multiply(B, disp);
			uB(0, 0) -= gpTE(0, 0);
			uB(1, 0) -= gpTE(0, 1);
			uB(2, 0) -= gpTE(0, 2);

			Sigma.multiply(C, uB, 1, 1);
		}
	}

	iterator.componentStress.data[0] = Sigma(0, 0) / N.size();
	iterator.componentStress.data[1] = Sigma(1, 0) / N.size();
	iterator.componentStress.data[2] = Sigma(2, 0) / N.size();
	iterator.componentStress.data[3] = Sigma(3, 0) / N.size();
	iterator.componentStress.data[4] = Sigma(4, 0) / N.size();
	iterator.componentStress.data[5] = Sigma(5, 0) / N.size();

	MATH::upDense3x3EigenValues(iterator.componentStress.data, iterator.principalStress.data);

	iterator.vonMisesStress.data[0] = sqrt((
			(iterator.principalStress.data[0] - iterator.principalStress.data[1]) *
			(iterator.principalStress.data[0] - iterator.principalStress.data[1]) +
			(iterator.principalStress.data[0] - iterator.principalStress.data[2]) *
			(iterator.principalStress.data[0] - iterator.principalStress.data[2]) +
			(iterator.principalStress.data[1] - iterator.principalStress.data[2]) *
			(iterator.principalStress.data[1] - iterator.principalStress.data[2])) / 2);
}

void StructuralMechanics3DKernel::nodeSolution(ElasticityElementIterator &iterator)
{

}

