
#include "structuralmechanics.h"
#include "assembler.hpp"

#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/boundaryregionstore.h"

#include "math/physics/matrix_distributed.h"

#include <numeric>
#include <algorithm>

#include <iostream>

namespace espreso {

NodeData* StructuralMechanics::Results::displacement = nullptr;
NodeData* StructuralMechanics::Results::thickness = nullptr;
ElementData* StructuralMechanics::Results::principalStress = nullptr;
ElementData* StructuralMechanics::Results::componentStress = nullptr;
ElementData* StructuralMechanics::Results::vonMisesStress = nullptr;
ElementData* StructuralMechanics::Results::isPlastized = nullptr;

StructuralMechanics::StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{
	elements.stiffness.setConstness(false);
	elements.stiffness.resize();
	elements.rhs.setConstness(false);
	elements.rhs.resize();
	elements.nrhs.setConstness(false);
	elements.nrhs.resize();
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		elements.boundary.rhs.regions[r].setConstness(false);
		elements.boundary.rhs.regions[r].resize();
	}

	subkernels.resize(info::mesh->elements->eintervals.size());
	boundary.resize(info::mesh->boundaryRegions.size());
	for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			boundary[r].resize(info::mesh->boundaryRegions[r]->eintervals.size());

		} else {
			boundary[r].resize(info::env::threads);
		}
	}

	for (int t = 0; t < info::env::threads; ++t) {
		for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
			for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
				subkernels[i].code = info::mesh->elements->eintervals[i].code;
				subkernels[i].elements = info::mesh->elements->eintervals[i].end - info::mesh->elements->eintervals[i].begin;
				subkernels[i].chunks = subkernels[i].elements / SIMD::size + (subkernels[i].elements % SIMD::size ? 1 : 0);
			}

			for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
				if (info::mesh->boundaryRegions[r]->dimension) {
					for (esint i = info::mesh->boundaryRegions[r]->eintervalsDistribution[d]; i < info::mesh->boundaryRegions[r]->eintervalsDistribution[d + 1]; ++i) {
						boundary[r][i].code = info::mesh->boundaryRegions[r]->eintervals[i].code;
						boundary[r][i].elements = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;
						boundary[r][i].chunks = boundary[r][i].elements / SIMD::size + (boundary[r][i].elements % SIMD::size ? 1 : 0);
					}
				} else {
					boundary[r][t].code = static_cast<int>(Element::CODE::POINT1);
					boundary[r][t].elements = info::mesh->boundaryRegions[r]->nodes->datatarray().size(t);
					boundary[r][t].chunks = boundary[r][t].elements / SIMD::size + (boundary[r][t].elements % SIMD::size ? 1 : 0);
				}
			}
		}
	}


	GaussPoints<Element::CODE::LINE2    ,  2, StructuralMechanicsGPC::LINE2    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE3,  3, StructuralMechanicsGPC::TRIANGLE3, 2>::set();
	GaussPoints<Element::CODE::SQUARE4  ,  4, StructuralMechanicsGPC::SQUARE4  , 2>::set();
	GaussPoints<Element::CODE::TETRA4   ,  4, StructuralMechanicsGPC::TETRA4   , 3>::set();
	GaussPoints<Element::CODE::PYRAMID5 ,  5, StructuralMechanicsGPC::PYRAMID5 , 3>::set();
	GaussPoints<Element::CODE::PRISMA6  ,  6, StructuralMechanicsGPC::PRISMA6  , 3>::set();
	GaussPoints<Element::CODE::HEXA8    ,  8, StructuralMechanicsGPC::HEXA8    , 3>::set();
	GaussPoints<Element::CODE::LINE3    ,  3, StructuralMechanicsGPC::LINE3    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE6,  6, StructuralMechanicsGPC::TRIANGLE6, 2>::set();
	GaussPoints<Element::CODE::SQUARE8  ,  8, StructuralMechanicsGPC::SQUARE8  , 2>::set();
	GaussPoints<Element::CODE::TETRA10  , 10, StructuralMechanicsGPC::TETRA10  , 3>::set();
	GaussPoints<Element::CODE::PYRAMID13, 13, StructuralMechanicsGPC::PYRAMID13, 3>::set();
	GaussPoints<Element::CODE::PRISMA15 , 15, StructuralMechanicsGPC::PRISMA15 , 3>::set();
	GaussPoints<Element::CODE::HEXA20   , 20, StructuralMechanicsGPC::HEXA20   , 3>::set();
}

void StructuralMechanics::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");

	validateRegionSettings("MATERIAL", settings.material_set);
	validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
	validateRegionSettings("THICKNESS", settings.thickness);

	if (Results::displacement == nullptr) {
		Results::displacement = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "DISPLACEMENT");
	}
	if (Results::thickness == nullptr && info::mesh->dimension == 2) {
		Results::thickness = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "THICKNESS");
	}
	if (info::ecf->output.results_selection.stress && Results::principalStress == nullptr) {
		Results::principalStress = info::mesh->elements->appendData(info::mesh->dimension    , NamedData::DataType::NUMBERED   , "PRINCIPAL_STRESS");
		Results::componentStress = info::mesh->elements->appendData(info::mesh->dimension * 2, NamedData::DataType::TENSOR_SYMM, "COMPONENT_STRESS");
		Results::vonMisesStress  = info::mesh->elements->appendData(                        1, NamedData::DataType::SCALAR     , "VON_MISES_STRESS");
	}

	for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
		if (info::mesh->materials[i]->material_model == MaterialConfiguration::MATERIAL_MODEL::PLASTICITY) {
			if (Results::isPlastized == nullptr) {
				Results::isPlastized = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "IS_PLASTIZED");
			}
		}
	}


	eslog::info(" ============================================================================================= \n");
	bool correct = true;
	if (configuration.displacement.size()) {
		correct &= checkBoundaryParameter("FIXED DISPLACEMENT", configuration.displacement);
//		generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 0, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][0][s] = value; });
//		generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 1, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][1][s] = value; });
//		generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 2, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][2][s] = value; });
	}

	if (step::step.loadstep == 0) {
		///////////////////////////////////// Set materials and check if there is not any incorrect region intersection
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		eslog::info("\n  MATERIALS                                                                                    \n");
		eslog::info(" --------------------------------------------------------------------------------------------- \n");
		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
			MaterialConfiguration *mat = info::mesh->materials[i];

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				eslog::info("    COORDINATE SYSTEM:                                                              CARTESIAN \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("ROTATION.Z", mat->coordinate_system.rotation.z);
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("ROTATION.X", mat->coordinate_system.rotation.x);
					correct &= checkExpression("ROTATION.Y", mat->coordinate_system.rotation.y);
					correct &= checkExpression("ROTATION.Z", mat->coordinate_system.rotation.z);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("SPHERICAL coordinate system is not supported in 2D.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
					correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
					correct &= checkExpression("CENTER.Z", mat->coordinate_system.center.z);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
					correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("CENTER.X", mat->coordinate_system.center.x);
					correct &= checkExpression("CENTER.Y", mat->coordinate_system.center.y);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			if (info::mesh->dimension == 2) {
				switch (settings.element_behaviour) {
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
					eslog::info("     ELEMENT BEHAVIOR:                                                           PLANE STRAIN \n");
					break;
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
					eslog::info("     ELEMENT BEHAVIOR:                                                           PLANE STRESS \n");
					break;
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
					eslog::info("     ELEMENT BEHAVIOR:                                            PLANE STRESS WITH THICKNESS \n");
					correct &= checkElementParameter("THICKNESS", settings.thickness);
					break;
				case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
					eslog::info("     ELEMENT BEHAVIOR:                                                           AXISYMMETRIC \n");
					correct &= checkElementParameter("THICKNESS", settings.thickness);
					break;
				}
				eslog::info("                                                                                               \n");
			}

			correct &= checkExpression("DENSITY", mat->density);
			correct &= checkExpression("HEAT CAPACITY", mat->heat_capacity);
			eslog::info("                                                                                               \n");

			switch (mat->material_model) {
			case MaterialConfiguration::MATERIAL_MODEL::PLASTICITY:
				eslog::info("     PLASTICITY MODEL:                                                               ISOTROPIC \n");
				eslog::info("                                                                                               \n");
				correct &= checkExpression("INITIAL_YIELD_STRESS", mat->plasticity_properties.initial_yield_stress);
				correct &= checkExpression("ISOTROPIC_HARDENING", mat->plasticity_properties.isotropic_hardening);
				correct &= checkExpression("KINEMATIC_HARDENING", mat->plasticity_properties.kinematic_hardening);
				eslog::info("                                                                                               \n");
				/* no break */
			case MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC:
				switch (mat->linear_elastic_properties.model) {
				case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
					eslog::info(" LINEAR ELASTIC MODEL:                                                              ISOTROPIC \n");
					correct &= checkExpression("EX", mat->linear_elastic_properties.young_modulus.get(0, 0));
					correct &= checkExpression("MI", mat->linear_elastic_properties.poisson_ratio.get(0, 0));
					break;
				case LinearElasticPropertiesConfiguration::MODEL::ORTHOTROPIC:
					eslog::info("                MODEL:                                                            ORTHOTROPIC \n");
					correct &= checkExpression("EX", mat->linear_elastic_properties.young_modulus.get(0, 0));
					correct &= checkExpression("EY", mat->linear_elastic_properties.young_modulus.get(1, 1));
					correct &= checkExpression("EZ", mat->linear_elastic_properties.young_modulus.get(2, 2));
					correct &= checkExpression("MIXY", mat->linear_elastic_properties.poisson_ratio.get(0, 0));
					correct &= checkExpression("MIXZ", mat->linear_elastic_properties.poisson_ratio.get(1, 1));
					correct &= checkExpression("MIYZ", mat->linear_elastic_properties.poisson_ratio.get(2, 2));
					correct &= checkExpression("GXY", mat->linear_elastic_properties.shear_modulus.get(0, 0));
					correct &= checkExpression("GXZ", mat->linear_elastic_properties.shear_modulus.get(1, 1));
					correct &= checkExpression("GYZ", mat->linear_elastic_properties.shear_modulus.get(2, 2));
					break;
				case LinearElasticPropertiesConfiguration::MODEL::ANISOTROPIC:
					eslog::info("                MODEL:                                                            ANISOTROPIC \n");
					break;
				}
				break;
			case MaterialConfiguration::MATERIAL_MODEL::HYPER_ELASTIC:
				eslog::info("                                                                                HYPER ELASTIC \n");
				break;
			}
			eslog::info("                                                                                               \n");
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		printMaterials(settings.material_set);
		eslog::info(" ============================================================================================= \n");
	}

	if (configuration.acceleration.size()) {
		correct &= checkElementParameter("ACCELERATION", configuration.acceleration);
	}
	if (configuration.angular_velocity.size()) {
		switch (info::mesh->dimension) {
		case 2:
			switch (settings.element_behaviour) {
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
				correct &= checkElementParameter("ANGULAR_VELOCITY.Z", configuration.angular_velocity, 2); break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
				correct &= checkElementParameter("ANGULAR_VELOCITY.Y", configuration.angular_velocity, 1); break;
			} break;
		case 3:
			correct &= checkElementParameter("ANGULAR_VELOCITY", configuration.angular_velocity); break;
		}
	}

	if (configuration.normal_pressure.size()) {
		correct &= checkBoundaryParameter("NORMAL_PRESSURE", configuration.normal_pressure);
	}

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
		bool rotated = mat->coordinate_system.isRotated() && mat->linear_elastic_properties.model != LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC;
		bool cartesian = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
		bool gpcoo = mat->linear_elastic_properties.needCoordinates() || getExpression(i, configuration.angular_velocity);
		gpcoo |= settings.element_behaviour == StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC;
		gpcoo |= mat->coordinate_system.type != CoordinateSystemConfiguration::TYPE::CARTESIAN;
//		bool gptemp = mat->linear_elastic_properties.needTemperature();
		esint eoffset = info::mesh->elements->eintervals[i].begin;

		if (info::mesh->dimension == 2) {
			subkernels[i].thickness.activate(getExpression(i, settings.thickness), info::mesh->elements->nodes->cbegin() + eoffset, info::mesh->elements->nodes->cend(), Results::thickness->data.data());
		}

		if (Results::principalStress) {
			subkernels[i].stress.activate(i, Results::principalStress, Results::componentStress, Results::vonMisesStress);
			subkernels[i].displacement.activate(info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin, info::mesh->elements->nodes->cend(), Results::displacement->data.data(), false);
			subkernels[i].elasticity.action |= Assembler::SOLUTION;
		}

		subkernels[i].coordinates.activate(info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin, info::mesh->elements->nodes->cend(), !cartesian || gpcoo);
		subkernels[i].elasticity.activate(settings.element_behaviour, &mat->linear_elastic_properties, rotated);
		if (mat->material_model == MaterialBaseConfiguration::MATERIAL_MODEL::PLASTICITY) {
			subkernels[i].elasticity.indirect = true;
			subkernels[i].plasticity.activate(i, settings.element_behaviour, &mat->plasticity_properties, Results::isPlastized, (elements.nrhs.data->begin() + i)->data());
			subkernels[i].displacement.activate(info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin, info::mesh->elements->nodes->cend(), Results::displacement->data.data(), true);
		}
		subkernels[i].coosystem.activate(mat->coordinate_system, subkernels[i].elasticity.isconst, rotated);
		subkernels[i].material.activate(mat);

		subkernels[i].K.activate((elements.stiffness.data->begin() + i)->data());
		subkernels[i].acceleration.activate(getExpression(i, configuration.acceleration), (elements.rhs.data->begin() + i)->data());
		subkernels[i].angularVelocity.activate(getExpression(i, configuration.angular_velocity), (elements.rhs.data->begin() + i)->data());
		if (Results::principalStress) {
			subkernels[i].stress.activate(i, Results::principalStress, Results::componentStress, Results::vonMisesStress);
		}
	}

	for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (info::mesh->boundaryRegions[r]->dimension) {
			for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				boundary[r][i].coordinates.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cend(), settings.element_behaviour == StructuralMechanicsGlobalSettings::ELEMENT_BEHAVIOUR::AXISYMMETRIC);
				boundary[r][i].normalPressure.activate(getExpression(info::mesh->boundaryRegions[r]->name, configuration.normal_pressure), (elements.boundary.rhs.regions[r].data->begin() + i)->data());
				boundary[r][i].integration.withNormal = boundary[r][i].normalPressure.isactive;
			}
		} else {
			for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
				boundary[r][t].coordinates.activate(region->nodes->cbegin(t), region->nodes->cend(), false);
			}
		}
	}

	for (auto it = configuration.displacement.begin(); it != configuration.displacement.end(); ++it) {
		size_t r = info::mesh->bregionIndex(it->first);
		for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
			boundary[r][t].displacement.activate(it->second);
		}
	}

	assemble(Action::PREPROCESS);
	size_t esize = 0;
	std::vector<double> volume(subkernels.size()), surface(boundary.size());
	for (size_t i = 0; i < subkernels.size(); ++i) {
		esize = std::max(subkernels[i].esize, esize);
		volume[i] = subkernels[i].volume;
	}
	for (size_t r = 1; r < boundary.size(); ++r) {
		for (size_t i = 0; i < boundary[r].size(); ++i) {
			surface[r] += boundary[r][i].surface;
		}
		for (size_t i = 0; i < boundary[r].size(); ++i) {
			boundary[r][i].surface = surface[r];
		}
	}
	printElementVolume(volume);
	printBoundarySurface(surface);

	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	eslog::info("  SIMD SIZE                                                                                 %lu \n", SIMD::size);
	eslog::info("  MAX ELEMENT SIZE                                                                   %6lu B \n", esize);
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
	} else {
		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
	}
	eslog::info(" ============================================================================================= \n");
}

void StructuralMechanics::connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		subkernels[i].Kfiller.activate(i, info::mesh->dimension, subkernels[i].elements, (elements.stiffness.data->begin() + i)->data(), K);
		subkernels[i].K.shape = K->shape;

		subkernels[i].RHSfiller.activate(i, info::mesh->dimension, subkernels[i].elements, (elements.rhs.data->begin() + i)->data(), f);
		subkernels[i].nRHSfiller.activate(i, info::mesh->dimension, subkernels[i].elements, (elements.nrhs.data->begin() + i)->data(), nf);
	}

	for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				boundary[r][i].RHSfiller.activate(r, i, info::mesh->dimension, boundary[r][i].elements, (elements.boundary.rhs.regions[r].data->begin() + i)->data(), f);
			}
		}
	}
	for (auto it = configuration.displacement.begin(); it != configuration.displacement.end(); ++it) {
		size_t r = info::mesh->bregionIndex(it->first);
		for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
			boundary[r][t].dirichlet.activate(r, t, info::mesh->dimension, boundary[r][t].elements, nullptr, dirichlet);
		}
	}
}

void StructuralMechanics::evaluate(const step::Step &step, const step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Matrix_Base<double> *C, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
	for (size_t i = 0; i < subkernels.size(); ++i) {
		for (size_t e = 0; e < subkernels[i].expressions.size(); ++e) {
			#pragma omp parallel for
			for (int t = 0; t < info::env::threads; ++t) {
				subkernels[i].expressions[e]->evaluator->getSubstep(t) = (step.substep + 1) / (double)step.substeps;
				subkernels[i].expressions[e]->evaluator->getTime(t) = time.current;
			}
		}
	}
	for (size_t i = 0; i < boundary.size(); ++i) {
		for (size_t j = 0; j < boundary[i].size(); ++j) {
			for (size_t e = 0; e < boundary[i][j].expressions.size(); ++e) {
				#pragma omp parallel for
				for (int t = 0; t < info::env::threads; ++t) {
					boundary[i][j].expressions[e]->evaluator->getSubstep(t) = (step.substep + 1) / (double)step.substeps;
					boundary[i][j].expressions[e]->evaluator->getTime(t) = time.current;
				}
			}
		}
	}
	reset(K, M, C, f, nf, dirichlet);
	assemble(Action::ASSEMBLE);
	assemble(Action::FILL);
	update(K, M, C, f, nf, dirichlet);
}

void StructuralMechanics::updateSolution(Vector_Base<double> *x)
{
	x->storeTo(Results::displacement->data);
	assemble(Action::SOLUTION);
}

void StructuralMechanics::nextIteration(Vector_Base<double> *x)
{
	x->storeTo(Results::displacement->data);
	assemble(Action::ITERATION);
}

void StructuralMechanics::run(Action action, size_t interval)
{
	switch (action) {
	case Action::PREPROCESS:
	case Action::FILL:
		runPreprocess(action, interval);
		break;
	default:
		switch (subkernels[interval].code) {
		case static_cast<size_t>(Element::CODE::TRIANGLE3): runGroup<Element::CODE::TRIANGLE3>(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::TRIANGLE6): runGroup<Element::CODE::TRIANGLE6>(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::SQUARE4  ): runGroup<Element::CODE::SQUARE4  >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::SQUARE8  ): runGroup<Element::CODE::SQUARE8  >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::TETRA4   ): runGroup<Element::CODE::TETRA4   >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::TETRA10  ): runGroup<Element::CODE::TETRA10  >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::PYRAMID5 ): runGroup<Element::CODE::PYRAMID5 >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::PYRAMID13): runGroup<Element::CODE::PYRAMID13>(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::PRISMA6  ): runGroup<Element::CODE::PRISMA6  >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::PRISMA15 ): runGroup<Element::CODE::PRISMA15 >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::HEXA8    ): runGroup<Element::CODE::HEXA8    >(action, interval, settings.element_behaviour); break;
		case static_cast<size_t>(Element::CODE::HEXA20   ): runGroup<Element::CODE::HEXA20   >(action, interval, settings.element_behaviour); break;
		}
	}
}

void StructuralMechanics::run(Action action, size_t region, size_t interval)
{
	runBoundary(action, region, interval);
}

}
