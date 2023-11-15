
#include "heattransfer.h"
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

NodeData* HeatTransfer::Results::temperature = nullptr;
NodeData* HeatTransfer::Results::initialTemperature = nullptr;
NodeData* HeatTransfer::Results::thickness = nullptr;
ElementData* HeatTransfer::Results::translationMotion = nullptr;
ElementData* HeatTransfer::Results::gradient = nullptr;
ElementData* HeatTransfer::Results::flux = nullptr;

HeatTransfer::HeatTransfer(HeatTransfer *previous, HeatTransferConfiguration &settings, HeatTransferLoadStepConfiguration &configuration)
: Assembler(settings), settings(settings), configuration(configuration)
{
	elements.stiffness.setConstness(false);
	elements.stiffness.resize();
	elements.rhs.setConstness(false);
	elements.rhs.resize();
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
						size_t elements = info::mesh->boundaryRegions[r]->eintervals[i].end - info::mesh->boundaryRegions[r]->eintervals[i].begin;
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


	GaussPoints<Element::CODE::LINE2    ,  2, HeatTransferGPC::LINE2    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE3,  3, HeatTransferGPC::TRIANGLE3, 2>::set();
	GaussPoints<Element::CODE::SQUARE4  ,  4, HeatTransferGPC::SQUARE4  , 2>::set();
	GaussPoints<Element::CODE::TETRA4   ,  4, HeatTransferGPC::TETRA4   , 3>::set();
	GaussPoints<Element::CODE::PYRAMID5 ,  5, HeatTransferGPC::PYRAMID5 , 3>::set();
	GaussPoints<Element::CODE::PRISMA6  ,  6, HeatTransferGPC::PRISMA6  , 3>::set();
	GaussPoints<Element::CODE::HEXA8    ,  8, HeatTransferGPC::HEXA8    , 3>::set();
	GaussPoints<Element::CODE::LINE3    ,  3, HeatTransferGPC::LINE3    , 1>::set();
	GaussPoints<Element::CODE::TRIANGLE6,  6, HeatTransferGPC::TRIANGLE6, 2>::set();
	GaussPoints<Element::CODE::SQUARE8  ,  8, HeatTransferGPC::SQUARE8  , 2>::set();
	GaussPoints<Element::CODE::TETRA10  , 10, HeatTransferGPC::TETRA10  , 3>::set();
	GaussPoints<Element::CODE::PYRAMID13, 13, HeatTransferGPC::PYRAMID13, 3>::set();
	GaussPoints<Element::CODE::PRISMA15 , 15, HeatTransferGPC::PRISMA15 , 3>::set();
	GaussPoints<Element::CODE::HEXA20   , 20, HeatTransferGPC::HEXA20   , 3>::set();
}

void HeatTransfer::analyze()
{
	double start = eslog::time();
	eslog::info("\n ============================================================================================= \n");

	validateRegionSettings("MATERIAL", settings.material_set);
	validateRegionSettings("INITIAL TEMPERATURE", settings.initial_temperature);
	validateRegionSettings("THICKNESS", settings.thickness);

	if (Results::temperature == nullptr) {
		Results::temperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
	}
	if (Results::initialTemperature == nullptr) {
		Results::initialTemperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");
	}
	if (Results::thickness == nullptr && info::mesh->dimension == 2) {
		Results::thickness = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "THICKNESS");
	}
	if (info::ecf->output.results_selection.translation_motions && Results::translationMotion == nullptr) {
		Results::translationMotion = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "TRANSLATION_MOTION");
	}
	if (info::ecf->output.results_selection.gradient && Results::gradient == nullptr) {
		Results::gradient = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "GRADIENT");
	}
	if (info::ecf->output.results_selection.flux && Results::flux == nullptr) {
		Results::flux = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FLUX");
	}

	eslog::info(" ============================================================================================= \n");
	bool correct = true;

	if (configuration.temperature.size()) {
		correct &= checkBoundaryParameter("FIXED TEMPERATURE ON BOUNDARIES", configuration.temperature);
	}

	if (settings.initial_temperature.size()) {
		correct &= checkElementParameter("INITIAL TEMPERATURE", settings.initial_temperature);
	}

	if (step::step.loadstep == 0) {
		if (info::mesh->dimension == 2) {
			correct &= checkElementParameter("THICKNESS", settings.thickness);
		}
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

			correct &= checkExpression("DENSITY", mat->density);
			correct &= checkExpression("HEAT CAPACITY", mat->heat_capacity);
			eslog::info("                                                                                               \n");

		switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
				correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
				}
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkExpression("KYZ", mat->thermal_conductivity.values.get(1, 2));
					correct &= checkExpression("KXZ", mat->thermal_conductivity.values.get(0, 2));
				}
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
				if (info::mesh->dimension == 2) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(1, 0));
				}
				if (info::mesh->dimension == 3) {
					correct &= checkExpression("KXX", mat->thermal_conductivity.values.get(0, 0));
					correct &= checkExpression("KYY", mat->thermal_conductivity.values.get(1, 1));
					correct &= checkExpression("KZZ", mat->thermal_conductivity.values.get(2, 2));
					correct &= checkExpression("KXY", mat->thermal_conductivity.values.get(0, 1));
					correct &= checkExpression("KYZ", mat->thermal_conductivity.values.get(1, 2));
					correct &= checkExpression("KXZ", mat->thermal_conductivity.values.get(0, 2));
					correct &= checkExpression("KYX", mat->thermal_conductivity.values.get(1, 0));
					correct &= checkExpression("KZY", mat->thermal_conductivity.values.get(2, 1));
					correct &= checkExpression("KZX", mat->thermal_conductivity.values.get(2, 0));
				}
				break;
			}
			eslog::info("                                                                                               \n");
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		printMaterials(settings.material_set);
		eslog::info(" ============================================================================================= \n");
	}

	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		const MaterialConfiguration *mat = info::mesh->materials[info::mesh->elements->eintervals[i].material];
		bool indirect = configuration.translation_motions.size() && settings.sigma > 0;
		bool rotated = mat->coordinate_system.isRotated() && mat->thermal_conductivity.model != ThermalConductivityConfiguration::MODEL::ISOTROPIC;
		bool cartesian = mat->coordinate_system.type == CoordinateSystemConfiguration::TYPE::CARTESIAN;
		bool gpcoo = mat->thermal_conductivity.needCoordinates();
		bool gptemp = mat->thermal_conductivity.needTemperature();
		esint eoffset = info::mesh->elements->eintervals[i].begin;

		if (info::mesh->dimension == 2) {
			subkernels[i].thickness.activate(getExpression(i, settings.thickness), info::mesh->elements->nodes->cbegin() + eoffset, info::mesh->elements->nodes->cend(), Results::thickness->data.data());
		}
		subkernels[i].material.activate(mat);

		subkernels[i].coordinates.activate(info::mesh->elements->nodes->cbegin() + eoffset, info::mesh->elements->nodes->cend(), !cartesian || gpcoo);
		subkernels[i].conductivity.activate(&mat->thermal_conductivity, rotated || indirect);
		if (indirect || rotated) {
			subkernels[i].coosystem.activate(mat->coordinate_system, mat->coordinate_system.isConst(), rotated);
		}
		subkernels[i].heatSource.activate(getExpression(i, configuration.heat_source), (elements.rhs.data->begin() + i)->data());
		subkernels[i].advection.activate(getExpression(i, configuration.translation_motions), (elements.stiffness.data->begin() + i)->data(), settings.sigma);

		subkernels[i].K.activate((elements.stiffness.data->begin() + i)->data());

		if (Results::gradient != nullptr) {
			subkernels[i].gradient.activate(i, Results::gradient);
		}
		if (Results::flux != nullptr) {
			subkernels[i].flux.activate(i, Results::flux);
			subkernels[i].conductivity.action |= Action::SOLUTION;
			subkernels[i].coosystem.action |= Action::SOLUTION;
		}

		if (Results::gradient != nullptr || Results::flux != nullptr || gptemp) {
			subkernels[i].temperature.activate(info::mesh->elements->nodes->cbegin() + info::mesh->elements->eintervals[i].begin, info::mesh->elements->nodes->cend(), Results::temperature->data.data(), gptemp);
			if (!gptemp) {
				subkernels[i].temperature.action = Action::SOLUTION;
			}
		}
	}

	for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		const BoundaryRegionStore *region = info::mesh->boundaryRegions[r];
		if (info::mesh->boundaryRegions[r]->dimension) {
			for(size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				if (info::mesh->dimension == 2) {
					boundary[r][i].thickness.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cend(), Results::thickness->data.data());
				}
				boundary[r][i].coordinates.activate(region->elements->cbegin() + region->eintervals[i].begin, region->elements->cend(), false);
				boundary[r][i].heatFlow.activate(getExpression(info::mesh->boundaryRegions[r]->name, configuration.heat_flow), nullptr);
				boundary[r][i].heatFlux.activate(getExpression(info::mesh->boundaryRegions[r]->name, configuration.heat_flux), nullptr);

				auto convection = configuration.convection.find(info::mesh->boundaryRegions[r]->name);
				if (convection != configuration.convection.end()) {
					boundary[r][i].htc.activate(&convection->second.heat_transfer_coefficient, nullptr);
					boundary[r][i].externalTemperature.activate(&convection->second.external_temperature, nullptr);
				}
				if (boundary[r][i].heatFlow.isactive | boundary[r][i].heatFlux.isactive | boundary[r][i].htc.isactive) {
					boundary[r][i].externalHeat.activate(boundary[r][i].surface, (elements.boundary.rhs.regions[r].data->begin() + i)->data());
				}
			}
		} else {
			for(size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
				boundary[r][t].coordinates.activate(region->nodes->cbegin(t), region->nodes->cend(), false);
			}
		}
	}

	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
		size_t r = info::mesh->bregionIndex(it->first);
		for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
			boundary[r][t].temperature.activate(it->second);
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
	}
	printElementVolume(volume);
	printBoundarySurface(surface);
	for (size_t r = 1; r < boundary.size(); ++r) {
		for (size_t i = 0; i < boundary[r].size(); ++i) {
			boundary[r][i].surface = boundary[r][i].externalHeat.area = surface[r];
		}
	}

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

void HeatTransfer::connect(Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
	for(size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		subkernels[i].Kfiller.activate(i, 1, subkernels[i].elements, (elements.stiffness.data->begin() + i)->data(), K);
		subkernels[i].Mfiller.activate(i, 1, subkernels[i].elements, (elements.stiffness.data->begin() + i)->data(), M);
		subkernels[i].RHSfiller.activate(i, 1, subkernels[i].elements, (elements.rhs.data->begin() + i)->data(), f);
		subkernels[i].nRHSfiller.activate(i, 1, subkernels[i].elements, (elements.rhs.data->begin() + i)->data(), nf);
		subkernels[i].K.shape = K->shape;
	}

	for(size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
		if (info::mesh->boundaryRegions[r]->dimension) {
			for (size_t i = 0; i < info::mesh->boundaryRegions[r]->eintervals.size(); ++i) {
				boundary[r][i].RHSfiller.activate(r, i, 1, boundary[r][i].elements, (elements.boundary.rhs.regions[r].data->begin() + i)->data(), f);
			}
		}
	}
	for (auto it = configuration.temperature.begin(); it != configuration.temperature.end(); ++it) {
		size_t r = info::mesh->bregionIndex(it->first);
		for (size_t t = 0; t < info::mesh->boundaryRegions[r]->nodes->threads(); ++t) {
			boundary[r][t].dirichlet.activate(r, t, 1, boundary[r][t].elements, nullptr, dirichlet);
		}
	}
}

void HeatTransfer::run(Action action, size_t interval)
{
	switch (action) {
	case Action::PREPROCESS:
	case Action::FILL:
		runPreprocess(action, interval);
		break;
	default:
		switch (info::mesh->dimension) {
		case 3: runVolume(action, interval); break;
		case 2: runPlane(action, interval); break;
		}
	}
}

void HeatTransfer::run(Action action, size_t region, size_t interval)
{
	runBoundary(action, region, interval);
}

void HeatTransfer::evaluate(step::Time &time, Matrix_Base<double> *K, Matrix_Base<double> *M, Vector_Base<double> *f, Vector_Base<double> *nf, Vector_Base<double> *dirichlet)
{
	reset(K, f, dirichlet);
	assemble(Action::ASSEMBLE);
	assemble(Action::FILL);
	update(K, f);
}

void HeatTransfer::updateSolution(Vector_Base<double> *x)
{
	x->storeTo(Results::temperature->data);
	assemble(Action::SOLUTION);
}

}

