
#include "structuralmechanics.h"
#include "structuralmechanics.generator.h"
#include "assembler.hpp"

#include "analysis/assembler/operators/info.h"
#include "analysis/assembler/operators/basis.h"
#include "analysis/assembler/operators/coordinates.h"
#include "analysis/assembler/operators/displacement.h"
#include "analysis/assembler/operators/expression.h"
#include "analysis/assembler/operators/integration.h"
#include "analysis/assembler/operators/structuralmechanics.f.h"
#include "analysis/assembler/operators/structuralmechanics.K.h"
#include "analysis/assembler/operators/filler.h"
#include "analysis/assembler/operators/stress.h"

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

#include <iostream>

namespace espreso {

NodeData* StructuralMechanics::Results::displacement = nullptr;
ElementData* StructuralMechanics::Results::principalStress = nullptr;
ElementData* StructuralMechanics::Results::componentStress = nullptr;
ElementData* StructuralMechanics::Results::vonMisesStress = nullptr;

StructuralMechanics::StructuralMechanics(StructuralMechanics *previous, StructuralMechanicsConfiguration &settings, StructuralMechanicsLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{
	axisymmetric = settings.element_behaviour == StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC;

	elements.stiffness.setConstness(false);
	elements.stiffness.resize();
	elements.rhs.setConstness(false);
	elements.rhs.resize();
	for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		elements.boundary.rhs.regions[r].setConstness(false);
		elements.boundary.rhs.regions[r].resize();
		if (configuration.normal_pressure.end() != configuration.normal_pressure.find(info::mesh->boundaryRegions[r]->name)) {
			bfilter[r] = 1;
		}
	}
	cossin_conditions.resize(info::mesh->elements->eintervals.size());
}

void StructuralMechanics::initParameters()
{
	if (Results::displacement == nullptr) {
		Results::displacement = info::mesh->nodes->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "DISPLACEMENT");
	}
	if (info::ecf->output.results_selection.stress && Results::principalStress == nullptr) {
		Results::principalStress = info::mesh->elements->appendData(info::mesh->dimension    , NamedData::DataType::NUMBERED   , "PRINCIPAL_STRESS");
		Results::componentStress = info::mesh->elements->appendData(info::mesh->dimension * 2, NamedData::DataType::TENSOR_SYMM, "COMPONENT_STRESS");
		Results::vonMisesStress  = info::mesh->elements->appendData(                        1, NamedData::DataType::SCALAR     , "VON_MISES_STRESS");
	}
}

bool StructuralMechanics::initDisplacement()
{
	bool correct = true;

	if (configuration.displacement.size()) {
		correct &= checkBoundaryParameter("FIXED DISPLACEMENT", configuration.displacement);
		generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 0, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][0][s] = value; });
		generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 1, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][1][s] = value; });
		generateBoundaryExpression<ExternalNodeExpression>(axisymmetric, boundaryOps, configuration.displacement, 2, [] (auto &element, const size_t &n, const size_t &s, const double &value) { element.displacement[n][2][s] = value; });
	}
	return correct;
}

void StructuralMechanics::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");

	validateRegionSettings("MATERIAL", settings.material_set);
	validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
	validateRegionSettings("THICKNESS", settings.thickness);

	initParameters();

	eslog::info(" ============================================================================================= \n");
	bool correct = true;
	correct &= initDisplacement();

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
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
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

			switch (mat->material_model) {
			case MaterialConfiguration::MATERIAL_MODEL::LINEAR_ELASTIC:
				eslog::info("                                                                               LINEAR ELASTIC \n");
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
				}
				eslog::info("                                                                                               \n");

				correct &= checkExpression("DENSITY", mat->density);
				correct &= checkExpression("HEAT CAPACITY", mat->heat_capacity);
				eslog::info("                                                                                               \n");

				switch (mat->linear_elastic_properties.model) {
				case LinearElasticPropertiesConfiguration::MODEL::ISOTROPIC:
					eslog::info("                MODEL:                                                              ISOTROPIC \n");
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

		switch (info::mesh->dimension) {
		case 2:
			if (settings.element_behaviour == StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC) {
				std::fill(etype.begin(), etype.end(), StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC);
			} else {
				std::fill(etype.begin(), etype.end(), StructuralMechanicsElementType::SYMMETRIC_PLANE);
			}
			break;
		case 3:
			std::fill(etype.begin(), etype.end(), StructuralMechanicsElementType::SYMMETRIC_VOLUME);
			break;
		}

		for (size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
			switch (info::mesh->boundaryRegions[r]->dimension) {
			case 0: std::fill(btype[r].begin(), btype[r].end(), StructuralMechanicsElementType::NODE); break;
			case 1: std::fill(btype[r].begin(), btype[r].end(), StructuralMechanicsElementType::EDGE); break;
			case 2: std::fill(btype[r].begin(), btype[r].end(), StructuralMechanicsElementType::FACE); break;
			}
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		printMaterials(settings.material_set);
		eslog::info(" ============================================================================================= \n");
	}

	generateBaseFunctions(etype, elementOps);
	generateBaseFunctions(axisymmetric, bfilter, boundaryOps);

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		auto procNodes = info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin;
		bool cooToGPs = settings.element_behaviour == StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC;

		if (getEvaluator(i, configuration.angular_velocity, 0) || getEvaluator(i, configuration.angular_velocity, 1) || getEvaluator(i, configuration.angular_velocity, 2)) {
			cooToGPs = true;
		}

		if (cooToGPs) {
			elementOps[i].push_back(generateElementOperator<CoordinatesToElementNodesAndGPs>(i, etype[i], procNodes));
		} else {
			elementOps[i].push_back(generateElementOperator<CoordinatesToElementNodes>(i, etype[i], procNodes));
		}

		if (Results::principalStress != nullptr) {
			elementOps[i].push_back(generateElementOperator<DisplacementToElementNodes>(i, etype[i], procNodes, Results::displacement->data.data()));
		}
	}

	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			bool cooToGPs = settings.element_behaviour == StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC;
			if (info::mesh->boundaryRegions[r]->dimension) {
				for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
					auto procNodes = info::mesh->boundaryRegions[r]->elements->cbegin() + info::mesh->boundaryRegions[r]->eintervals[i].begin;
					if (cooToGPs) {
						boundaryOps[r][i].push_back(generateBoundaryOperator<CoordinatesToElementNodesAndGPs>(axisymmetric, r, i, procNodes));
					} else {
						boundaryOps[r][i].push_back(generateBoundaryOperator<CoordinatesToElementNodes>(axisymmetric, r, i, procNodes));
					}
				}
			}
		}
	}

	if (info::mesh->dimension == 2) {
		generateElementExpression2D<ExternalGpsExpressionWithCoordinates>(etype, elementOps, settings.thickness, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.thickness[gp][s] = value; });
	}

	generateElasticity();
	generateElementOperators<Integration>(etype, elementOps);
	if (configuration.normal_pressure.size()) {
		generateBoundaryOperators<IntegrationWithNormal>(axisymmetric, bfilter, boundaryOps);
	} else {
		generateBoundaryOperators<Integration>(axisymmetric, bfilter, boundaryOps);
	}
	volume();

	generateElementOperators<StructuralMechanicsStiffness>(etype, elementOps, elements.stiffness);

	if (configuration.acceleration.size()) {
		correct &= checkElementParameter("ACCELERATION", configuration.acceleration);
		generateElementExpression<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.acceleration, 0, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.acceleration[gp][0][s] = value; });
		generateElementExpression<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.acceleration, 1, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.acceleration[gp][1][s] = value; });
		if (info::mesh->dimension == 3) {
			generateElementExpression<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.acceleration, 2, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.acceleration[gp][2][s] = value; });
		}
		generateElementOperators<Acceleration>(etype, elementOps, elements.rhs);
	}

	if (configuration.angular_velocity.size()) {
		switch (info::mesh->dimension) {
		case 3:
			correct &= checkElementParameter("ANGULAR_VELOCITY", configuration.angular_velocity);
			generateElementExpression3D<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.angular_velocity, 0, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angularVelocity[gp][0][s] = value; });
			generateElementExpression3D<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.angular_velocity, 1, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angularVelocity[gp][1][s] = value; });
			generateElementExpression3D<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.angular_velocity, 2, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angularVelocity[gp][2][s] = value; });
			break;
		case 2:
			switch (settings.element_behaviour) {
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRAIN:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS:
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::PLANE_STRESS_WITH_THICKNESS:
				correct &= checkElementParameter("ANGULAR_VELOCITY.Z", configuration.angular_velocity, 2);
				generateElementExpression2D<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.angular_velocity, 2, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angularVelocity[gp][s] = value; });
				break;
			case StructuralMechanicsConfiguration::ELEMENT_BEHAVIOUR::AXISYMMETRIC:
				correct &= checkElementParameter("ANGULAR_VELOCITY.Y", configuration.angular_velocity, 1);
				generateElementExpression2D<ExternalGpsExpressionWithCoordinates>(etype, elementOps, configuration.angular_velocity, 1, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.angularVelocity[gp][s] = value; });
				break;
			}
		}
		generateElementOperators<AngularVelocity>(etype, elementOps, elements.rhs);
	}

	if (configuration.normal_pressure.size()) {
		correct &= checkBoundaryParameter("NORMAL PRESSURE", configuration.normal_pressure);
		generateBoundaryExpression<ExternalGpsExpressionWithCoordinates>(axisymmetric, boundaryOps, configuration.normal_pressure, [] (auto &element, const size_t &gp, const size_t &s, const double &value) { element.ecf.normalPressure[gp][s] = value; });
	}
	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (bfilter[r]) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				boundaryOps[r][i].push_back(generateBoundaryOperator<NormalPressure>(axisymmetric, r, i, elements.boundary.rhs.regions[r]));
			}
		}
	}

	if (Results::principalStress != nullptr) {
		generateElementOperators<Stress>(etype, elementOps, Results::principalStress, Results::componentStress, Results::vonMisesStress);
	}

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		if (Results::principalStress == nullptr) {
			for (size_t j = 0; j < elementOps[i].size(); ++j) {
				ActionOperator::removeSolution(elementOps[i][j]->action);
			}
		}
	}


	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	eslog::info("  SIMD SIZE                                                                                 %lu \n", SIMD::size);
	eslog::info("  MAX ELEMENT SIZE                                                                   %6lu B \n", esize());
	eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	if (correct) {
		eslog::info("  PHYSICS CONFIGURED                                                               %8.3f s \n", eslog::time() - start);
	} else {
		eslog::globalerror("  PHYSICS CONFIGURATION FAILED                                                         \n");
	}
	eslog::info(" ============================================================================================= \n");
}

void StructuralMechanics::connect(SteadyState &scheme)
{
	this->K = scheme.K;
	this->f = scheme.f;
	switch (scheme.K->shape) {
	case Matrix_Shape::FULL:
		for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
			switch (etype[i]) {
			case StructuralMechanicsElementType::SYMMETRIC_PLANE:
			case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC:
			case StructuralMechanicsElementType::SYMMETRIC_VOLUME:
				elementOps[i].push_back(generateElementOperator<SymmetricToFullMatricFiller>(i, etype[i], info::mesh->dimension, elements.stiffness, scheme.K)); break;
				break;
			default:
				elementOps[i].push_back(generateElementOperator<GeneralMatricFiller>(i, etype[i], info::mesh->dimension, elements.stiffness, scheme.K)); break;
				break;
			}
		}
		break;
	case Matrix_Shape::UPPER: generateElementOperators<SymmetricMatricFiller>(etype, elementOps, info::mesh->dimension, elements.stiffness, scheme.K); break;
	}

	generateElementOperators<VectorFiller>(etype, elementOps, info::mesh->dimension, elements.rhs, scheme.f);

	for(size_t r = 0; r < info::mesh->boundaryRegions.size(); ++r) {
		if (configuration.displacement.end() == configuration.displacement.find(info::mesh->boundaryRegions[r]->name)) {
			if (bfilter[r]) {
				switch (info::mesh->boundaryRegions[r]->dimension) {
				case 0:
					for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
						boundaryOps[r][t].push_back(generateNodeFiller<VectorFiller>(r, t, info::mesh->dimension, elements.rhs, scheme.f));
					}
					break;
				case 1:
					for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
						boundaryOps[r][i].push_back(generateEdgeFiller<VectorFiller>(axisymmetric, r, i, info::mesh->dimension, elements.boundary.rhs.regions[r], scheme.f));
					}
					break;
				case 2:
					for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
						boundaryOps[r][i].push_back(generateFaceFiller<VectorFiller>(r, i, info::mesh->dimension, elements.boundary.rhs.regions[r], scheme.f));
					}
					break;
				}
			}
		} else {
			// DIRICHLET
			for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
				boundaryOps[r][t].push_back(generateNodeSetter<VectorSetter>(r, t, info::mesh->dimension, scheme.dirichlet, [] (auto &element, const size_t &n, const size_t &dof, const size_t &s) { return element.displacement[n][dof][s]; }));
			}
		}
	}
}

void StructuralMechanics::evaluate(SteadyState &scheme)
{
	reset(scheme.K, scheme.f, scheme.dirichlet);
	eslog::info("       = SIMD LOOP                                                      %12.8f s = \n", assemble(ActionOperator::Action::ASSEMBLE));
	eslog::info("       = FILL MATRICES                                                  %12.8f s = \n", assemble(ActionOperator::Action::FILL));
	update(scheme.K, scheme.f);
}

void StructuralMechanics::dryrun()
{
	if (this->K == nullptr) {
		this->K = new Matrix_Distributed<Matrix_CSR, double>();
		this->K->mapping.elements.resize(info::mesh->elements->eintervals.size());

		this->f = new Vector_Distributed<Vector_Dense, double>();
		this->f->mapping.elements.resize(info::mesh->elements->eintervals.size());
	}
	Assembler::measurements reassemble_time = Assembler::measurements();
	Assembler::measurements   assemble_time = Assembler::measurements();
	Assembler::measurements   solution_time = Assembler::measurements();
	int numreps = 10;
	for(int reps = 0; reps < numreps; reps++) {
		assemble_time   += assemble(ActionOperator::Action::ASSEMBLE);
		reassemble_time += assemble(ActionOperator::Action::REASSEMBLE);
		solution_time   += assemble(ActionOperator::Action::SOLUTION);
	}

	assemble_time.coreTime   /= static_cast<double>(numreps);
	reassemble_time.coreTime /= static_cast<double>(numreps);
	solution_time.coreTime   /= static_cast<double>(numreps);

	assemble_time.preprocessTime   /= static_cast<double>(numreps);
	reassemble_time.preprocessTime /= static_cast<double>(numreps);
	solution_time.preprocessTime   /= static_cast<double>(numreps);

	eslog::info("       = SIMD LOOP ASSEMBLE                                             %12.8f s = \n",   assemble_time.preprocessTime);
	std::cout<<"SCALING: "<<assemble_time.preprocessTime<<std::endl;
	eslog::info("       = SIMD LOOP ASSEMBLE                                             %12.8f s = \n",   assemble_time.coreTime);
	std::cout<<"SCALING: "<<assemble_time.coreTime<<std::endl;
	eslog::info("       = SIMD LOOP REASSEMBLE                                           %12.8f s = \n", reassemble_time.preprocessTime);
	std::cout<<"SCALING: "<<reassemble_time.preprocessTime<<std::endl;
	eslog::info("       = SIMD LOOP REASSEMBLE                                           %12.8f s = \n", reassemble_time.coreTime);
	std::cout<<"SCALING: "<<reassemble_time.coreTime<<std::endl;
	eslog::info("       = SIMD LOOP SOLUTION                                             %12.8f s = \n", solution_time.preprocessTime);
	std::cout<<"SCALING: "<<solution_time.preprocessTime<<std::endl;
	eslog::info("       = SIMD LOOP SOLUTION                                             %12.8f s = \n", solution_time.coreTime);
	std::cout<<"SCALING: "<<solution_time.coreTime<<std::endl;
}


void StructuralMechanics::volume()
{
	std::vector<double> evolume(info::mesh->elements->eintervals.size());
	std::vector<double> bvolume(info::mesh->boundaryRegions.size());

	generateElementOperators<Volume>(etype, elementOps, evolume);
	generateBoundaryOperators<Volume>(axisymmetric, bfilter, boundaryOps, bvolume);
	assemble(ActionOperator::Action::ASSEMBLE);
	dropLastOperators(elementOps);
	dropLastOperators(bfilter, boundaryOps);

	printElementVolume(evolume);
	printBoundarySurface(bvolume);
}

size_t StructuralMechanics::esize()
{
	size_t max = 0;
	std::vector<size_t> esize(info::mesh->elements->eintervals.size());

	generateElementOperators<DataDescriptorElementSize>(etype, elementOps, esize);
	for (size_t i = 0; i < elementOps.size(); ++i) {
		max = std::max(max, esize[i]);
	}
	dropLastOperators(elementOps);
	return max;
}

void StructuralMechanics::updateSolution(SteadyState &scheme)
{
	scheme.x->storeTo(Results::displacement->data);
	assemble(ActionOperator::Action::SOLUTION);
}

template <>
Assembler::measurements StructuralMechanics::instantiate2D<StructuralMechanicsElementType::NODE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	return loop<StructuralMechanicsDataDescriptor, 1, StructuralMechanicsGPC::POINT1, 2, 0, StructuralMechanicsElementType::NODE>(action, ops, elements);
}

template <>
Assembler::measurements StructuralMechanics::instantiate2D<StructuralMechanicsElementType::EDGE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2): return loop<StructuralMechanicsDataDescriptor, 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE>(action, ops,elements); break;
	case static_cast<size_t>(Element::CODE::LINE3): return loop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE>(action, ops,elements); break;
	default: return measurements();
	}
}

template <>
Assembler::measurements StructuralMechanics::instantiate2D<StructuralMechanicsElementType::EDGE_AXISYMMETRIC>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2): return loop<StructuralMechanicsDataDescriptor, 2, StructuralMechanicsGPC::LINE2, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3): return loop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::LINE3, 2, 1, StructuralMechanicsElementType::EDGE_AXISYMMETRIC>(action, ops, elements); break;
	default: return measurements();
	}
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiate2D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return loop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::TRIANGLE3, 2, 2, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return loop<StructuralMechanicsDataDescriptor, 6, StructuralMechanicsGPC::TRIANGLE6, 2, 2, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return loop<StructuralMechanicsDataDescriptor, 4, StructuralMechanicsGPC::SQUARE4  , 2, 2, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return loop<StructuralMechanicsDataDescriptor, 8, StructuralMechanicsGPC::SQUARE8  , 2, 2, etype>(action, ops, elements); break;
	default: return measurements();
	}
}

template <>
Assembler::measurements StructuralMechanics::instantiate3D<StructuralMechanicsElementType::NODE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	return loop<StructuralMechanicsDataDescriptor, 1, StructuralMechanicsGPC::POINT1, 3, 0, StructuralMechanicsElementType::NODE>(action, ops, elements);
}

template <>
Assembler::measurements StructuralMechanics::instantiate3D<StructuralMechanicsElementType::EDGE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::LINE2): return loop<StructuralMechanicsDataDescriptor, 2, StructuralMechanicsGPC::LINE2, 3, 1, StructuralMechanicsElementType::EDGE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::LINE3): return loop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::LINE3, 3, 1, StructuralMechanicsElementType::EDGE>(action, ops, elements); break;
	default: return measurements();
	};
}

template <>
Assembler::measurements StructuralMechanics::instantiate3D<StructuralMechanicsElementType::FACE>(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TRIANGLE3): return loop<StructuralMechanicsDataDescriptor, 3, StructuralMechanicsGPC::TRIANGLE3, 3, 2, StructuralMechanicsElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::TRIANGLE6): return loop<StructuralMechanicsDataDescriptor, 6, StructuralMechanicsGPC::TRIANGLE6, 3, 2, StructuralMechanicsElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE4):   return loop<StructuralMechanicsDataDescriptor, 4, StructuralMechanicsGPC::SQUARE4  , 3, 2, StructuralMechanicsElementType::FACE>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::SQUARE8):   return loop<StructuralMechanicsDataDescriptor, 8, StructuralMechanicsGPC::SQUARE8  , 3, 2, StructuralMechanicsElementType::FACE>(action, ops, elements); break;
	default: return measurements();
	}
}

template <int etype>
Assembler::measurements StructuralMechanics::instantiate3D(ActionOperator::Action action, int code, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (code) {
	case static_cast<size_t>(Element::CODE::TETRA4):    return loop<StructuralMechanicsDataDescriptor,  4, StructuralMechanicsGPC::TETRA4    , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::TETRA10):   return loop<StructuralMechanicsDataDescriptor, 10, StructuralMechanicsGPC::TETRA10   , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID5):  return loop<StructuralMechanicsDataDescriptor,  5, StructuralMechanicsGPC::PYRAMID5  , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PYRAMID13): return loop<StructuralMechanicsDataDescriptor, 13, StructuralMechanicsGPC::PYRAMID13 , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA6):   return loop<StructuralMechanicsDataDescriptor,  6, StructuralMechanicsGPC::PRISMA6   , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::PRISMA15):  return loop<StructuralMechanicsDataDescriptor, 15, StructuralMechanicsGPC::PRISMA15  , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA8):     return loop<StructuralMechanicsDataDescriptor,  8, StructuralMechanicsGPC::HEXA8     , 3, 3, etype>(action, ops, elements); break;
	case static_cast<size_t>(Element::CODE::HEXA20):    return loop<StructuralMechanicsDataDescriptor, 20, StructuralMechanicsGPC::HEXA20    , 3, 3, etype>(action, ops, elements); break;
	default: return measurements();
	}
}

Assembler::measurements StructuralMechanics::instantiate(ActionOperator::Action action, int code, int etype, const std::vector<ActionOperator*> &ops, size_t interval, esint elements)
{
	switch (info::mesh->dimension) {
	case 2:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_PLANE             : return instantiate2D<StructuralMechanicsElementType::SYMMETRIC_PLANE             >(action, code, ops, interval, elements);
		case StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC: return instantiate2D<StructuralMechanicsElementType::SYMMETRIC_PLANE_AXISYMMETRIC>(action, code, ops, interval, elements);

		// boundary
		case StructuralMechanicsElementType::EDGE:
			if (axisymmetric) {
				return instantiate2D<StructuralMechanicsElementType::EDGE_AXISYMMETRIC>(action, code, ops, interval, elements);
			} else {
				return instantiate2D<StructuralMechanicsElementType::EDGE>(action, code, ops, interval, elements);
			}
		case StructuralMechanicsElementType::NODE: return instantiate2D<StructuralMechanicsElementType::NODE>(action, code, ops, interval, elements);
		}
	case 3:
		switch (etype) {
		// elements
		case StructuralMechanicsElementType::SYMMETRIC_VOLUME: return instantiate3D<StructuralMechanicsElementType::SYMMETRIC_VOLUME>(action, code, ops, interval, elements);

		// boundary
		case StructuralMechanicsElementType::FACE: return instantiate3D<StructuralMechanicsElementType::FACE>(action, code, ops, interval, elements);
		case StructuralMechanicsElementType::EDGE: return instantiate3D<StructuralMechanicsElementType::EDGE>(action, code, ops, interval, elements);
		case StructuralMechanicsElementType::NODE: return instantiate3D<StructuralMechanicsElementType::NODE>(action, code, ops, interval, elements);
		}
	}
	return measurements();
}

}
