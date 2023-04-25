
#ifndef SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_SUBKERNELS_H_
#define SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_SUBKERNELS_H_

#include "analysis/assembler/subkernel/subkernel.h"
#include "config/ecf/material/linearelasticproperties.h"

namespace espreso {

struct StructuralMechanicsGPC {
	enum: int {
		POINT1    =  1,

		LINE2     =  2,
		LINE3     =  3,

		TRIANGLE3 =  6,
		SQUARE4   =  4,

		TRIANGLE6 =  6,
		SQUARE8   =  9,

		TETRA4    =  4,
		PYRAMID5  =  8,
		PRISMA6   =  9,
		HEXA8     =  8,

		TETRA10   = 15,
		PYRAMID13 = 14,
		PRISMA15  =  9,
		HEXA20    =  8,
	};
};

enum class Behaviour {
	PLANE,
	AXISYMMETRIC,
	VOLUME
};

enum class ElasticityModel {
	ISOTROPIC,
	ORTHOTROPIC,
	SYMMETRIC,
	ANISOTROPIC
};

template <size_t gps, size_t ndim, enum Behaviour, enum ElasticityModel> struct StructuralElasticity;

template <size_t gps, size_t ndim, enum Behaviour behaviour> struct StructuralElasticity<gps, ndim, behaviour, ElasticityModel::ISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][3];
};

template <size_t gps> struct StructuralElasticity<gps, 2, Behaviour::PLANE, ElasticityModel::ORTHOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][6]; // TODO
};

template <size_t gps> struct StructuralElasticity<gps, 2, Behaviour::PLANE, ElasticityModel::SYMMETRIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][6];
};

template <size_t gps> struct StructuralElasticity<gps, 2, Behaviour::PLANE, ElasticityModel::ANISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][9];
};

template <size_t gps> struct StructuralElasticity<gps, 2, Behaviour::AXISYMMETRIC, ElasticityModel::ORTHOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][10]; // TODO
};

template <size_t gps> struct StructuralElasticity<gps, 2, Behaviour::AXISYMMETRIC, ElasticityModel::SYMMETRIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][10];
};

template <size_t gps> struct StructuralElasticity<gps, 2, Behaviour::AXISYMMETRIC, ElasticityModel::ANISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][16];
};

template <size_t gps> struct StructuralElasticity<gps, 3, Behaviour::VOLUME, ElasticityModel::ORTHOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][9];
};

template <size_t gps> struct StructuralElasticity<gps, 3, Behaviour::VOLUME, ElasticityModel::SYMMETRIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][21];
};

template <size_t gps> struct StructuralElasticity<gps, 3, Behaviour::VOLUME, ElasticityModel::ANISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD elasticity[gps][36];
};

template <size_t gps, size_t ndim> struct ElasticityRotation {
	alignas(SIMD::size * sizeof(double)) SIMD center[gps][ndim];
};

template <size_t gps, size_t ndim> struct ElasticityRotationMatrix {
	alignas(SIMD::size * sizeof(double)) SIMD cossin[gps][4 * ndim];
};

template <size_t gps, size_t ndim> struct StructuralMechanicsElementParameters;

template <size_t gps> struct StructuralMechanicsElementParameters<gps, 2> {
	alignas(SIMD::size * sizeof(double)) SIMD thickness   [gps];

	alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];

	alignas(SIMD::size * sizeof(double)) SIMD acceleration   [gps][2];
	alignas(SIMD::size * sizeof(double)) SIMD angularVelocity[gps][1];

	StructuralMechanicsElementParameters()
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			thickness[gp] = load1(1);
		}
	}
};

template <size_t gps> struct StructuralMechanicsElementParameters<gps, 3> {
	alignas(SIMD::size * sizeof(double)) SIMD density     [gps];
	alignas(SIMD::size * sizeof(double)) SIMD heatCapacity[gps];

	alignas(SIMD::size * sizeof(double)) SIMD acceleration   [gps][3];
	alignas(SIMD::size * sizeof(double)) SIMD angularVelocity[gps][3];
};

template <size_t gps, size_t ndim, ElasticityModel model> struct ElasticityParameters {
	alignas(SIMD::size * sizeof(double)) SIMD youngModulus[gps][ndim];
	alignas(SIMD::size * sizeof(double)) SIMD poissonRatio[gps][ndim];
	alignas(SIMD::size * sizeof(double)) SIMD shearModulus[gps][ndim];
};

template <size_t gps, size_t ndim> struct ElasticityParameters<gps, ndim, ElasticityModel::ISOTROPIC> {
	alignas(SIMD::size * sizeof(double)) SIMD youngModulus[gps];
	alignas(SIMD::size * sizeof(double)) SIMD poissonRatio[gps];
};

template <size_t gps, size_t ndim> struct StructuralMechanicsBoundaryParameters;

template <size_t gps> struct StructuralMechanicsBoundaryParameters<gps, 2> {
	alignas(SIMD::size * sizeof(double)) SIMD thickness[gps];

	alignas(SIMD::size * sizeof(double)) SIMD normalPressure[gps];

	StructuralMechanicsBoundaryParameters()
	{
		for (size_t gp = 0; gp < gps; ++gp) {
			thickness[gp] = load1(1);
		}
	}
};

template <size_t gps> struct StructuralMechanicsBoundaryParameters<gps, 3> {
	alignas(SIMD::size * sizeof(double)) SIMD normalPressure[gps];
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim, enum Behaviour behaviour, enum ElasticityModel ecfmodel, enum ElasticityModel model> struct StructuralMechanicsElementDescriptor {
	virtual ~StructuralMechanicsElementDescriptor() {}

	struct Element:
			ElementCoordinates<nodes, gps, ndim>,
			ElementTemperature<nodes, gps>,
			ElementIntegration<nodes, gps, edim>,
			StructuralElasticity<gps, ndim, behaviour, model>,
			ElasticityRotationMatrix<gps, ndim>
	{
		struct ECF:
				StructuralElasticity<gps, ndim, behaviour, ecfmodel>,
				ElasticityRotation<gps, ndim>,
				StructuralMechanicsElementParameters<gps, ndim>,
				ElasticityParameters<gps, ndim, ecfmodel>
		{

		} ecf;
	};

	virtual void simd(Element &element) =0;
};

template <size_t nodes, size_t gps, size_t ndim, size_t edim> struct StructuralMechanicsBoundaryDescriptor {
	virtual ~StructuralMechanicsBoundaryDescriptor() {}

	struct Element:
			ElementCoordinates<nodes, gps, ndim>,
			ElementTemperature<nodes, gps>,
			ElementDisplacement<nodes, ndim>,
			ElementIntegration<nodes, gps, edim>,
			BondaryNormal<gps, ndim>
	{
		struct ECF:
				StructuralMechanicsBoundaryParameters<gps, ndim>
		{

		} ecf;
	};

	virtual void simd(Element &element) =0;
};

}



#endif /* SRC_ANALYSIS_ASSEMBLER_SUBKERNEL_STRUCTURALMECHANICS_SUBKERNELS_H_ */
