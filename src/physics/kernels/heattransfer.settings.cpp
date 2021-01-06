
#include "heattransfer.kernel.h"
#include "solverdataprovider/heattransfer.provider.h"

#include "basis/evaluator/evaluator.h"
#include "basis/utilities/utils.h"
#include "basis/utilities/parser.h"
#include "esinfo/meshinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "mesh/store/elementsregionstore.h"

#include "basis/utilities/print.h"

using namespace espreso;

NodeData* HeatTransferKernel::initialTemperature = NULL;
NodeData* HeatTransferKernel::temperature = NULL;
ElementData* HeatTransferKernel::translationMotion = NULL;
ElementData* HeatTransferKernel::phase = NULL;
ElementData* HeatTransferKernel::latentHeat = NULL;
ElementData* HeatTransferKernel::gradient = NULL;
ElementData* HeatTransferKernel::flux = NULL;

void HeatTransferKernel::createParameters()
{
	initialTemperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "INITIAL_TEMPERATURE");
	temperature = info::mesh->nodes->appendData(1, NamedData::DataType::SCALAR, "TEMPERATURE");
	if (info::ecf->output.results_selection.translation_motions) {
		translationMotion = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "TRANSLATION_MOTION");
	}
	if (info::ecf->output.results_selection.phase) {
		phase = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "PHASE");
	}
	if (info::ecf->output.results_selection.latent_heat) {
		latentHeat = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "LATENT_HEAT");
	}
	if (info::ecf->output.results_selection.gradient) {
		gradient = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "GRADIENT");
	}
	if (info::ecf->output.results_selection.flux) {
		flux = info::mesh->elements->appendData(info::mesh->dimension, NamedData::DataType::VECTOR, "FLUX");
	}
}

//ElementParameterInfo _thickness;
//	ElementParameterInfo _invJ, _detJ, _dND;
//	ElementParameterInfo _csCartesian, _csSpherical, _csCylindrical, _dens, _hc, _conductivity[4];
//	ElementParameterInfo _K, _isotropicK;
//	ElementParameterInfo _translationMotion, _u, _advection;

HeatTransferKernel::HeatTransferKernel(HeatTransferKernel *previous, HeatTransferLoadStepConfiguration &configuration)
: Kernel(new HeatTransferSolverDataProvider(configuration)),
  configuration(configuration),
  _N(1, ParameterInfo::Range::EACH_NODE_GP), _dN(info::mesh->dimension, ParameterInfo::Range::PER_NODE_GP), _w(1, ParameterInfo::Range::EACH_GP, 0),
  _time(1, ParameterInfo::Range::AGREGATED), _timestep(1, ParameterInfo::Range::AGREGATED),
  _ncoordinates(info::mesh->dimension, ParameterInfo::Range::PER_NODE),
  _ntemperature(1, ParameterInfo::Range::PER_NODE),
  _ninitTemperature(1, ParameterInfo::Range::PER_NODE, 273.15),
  _coordinates(info::mesh->dimension, ParameterInfo::Range::PER_GP),
  _initTemperature(1, ParameterInfo::Range::PER_GP), _temperature(1, ParameterInfo::Range::PER_GP),
  _thickness(1, ParameterInfo::Range::PER_GP, 1),
  _invJ(4, ParameterInfo::Range::PER_GP), _detJ(1, ParameterInfo::Range::PER_GP), _dND(info::mesh->dimension, ParameterInfo::Range::PER_NODE_GP),
  _csCartesian(info::mesh->dimension == 2 ? 1 : 3, ParameterInfo::Range::PER_GP, 0), _csSpherical(info::mesh->dimension, ParameterInfo::Range::PER_GP, 0), _csCylindrical(2, ParameterInfo::Range::PER_GP, 0),
  _dens(1, ParameterInfo::Range::PER_GP, 1), _hc(1, ParameterInfo::Range::PER_GP, 1),
  _conductivity{
	ElementParameterInfo{1, ParameterInfo::Range::PER_GP, 1},
	ElementParameterInfo{info::mesh->dimension, ParameterInfo::Range::PER_GP, 1},
	ElementParameterInfo{info::mesh->dimension + info::mesh->dimension / 2, ParameterInfo::Range::PER_GP, 1},
	ElementParameterInfo{info::mesh->dimension * info::mesh->dimension, ParameterInfo::Range::PER_GP, 1}
  },
  _K(4, ParameterInfo::Range::PER_GP, 0), _isotropicK(1, ParameterInfo::Range::PER_GP, 1),
  _translationMotion(info::mesh->dimension, ParameterInfo::Range::PER_GP, 0), _u(info::mesh->dimension, ParameterInfo::Range::PER_GP, 0), _advection(info::mesh->dimension, ParameterInfo::Range::PER_NODE_GP, 0),
  _stiffness(1, ParameterInfo::Range::PER_GP_GP), _mass(1, ParameterInfo::Range::PER_GP_GP),
  _btemperature(1, ParameterInfo::Range::PER_GP), _heatflow(1, ParameterInfo::Range::PER_GP), _heatflux(1, ParameterInfo::Range::PER_GP)
{
	solutions.push_back(VectorDense(temperature->data.size(), temperature->data.data()));

	_N.smartResize();
	_dN.smartResize();
	_w.smartResize();

	{
		int index = 0;
		for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei, ++index) {
			double *n = (_N.data->begin() + index)->data();
			double *dn = (_dN.data->begin() + index)->data();
			double *w = (_w.data->begin() + index)->data();
			esint nsize = (*Mesh::edata[ei->code].N)[0].ncols * (*Mesh::edata[ei->code].N)[0].nrows;
			for (size_t gp = 0; gp < Mesh::edata[ei->code].N->size(); ++gp) {
				memcpy(n + gp * nsize, (*Mesh::edata[ei->code].N)[gp].vals, sizeof(double) * nsize);
				memcpy(dn + info::mesh->dimension * gp * nsize, (*Mesh::edata[ei->code].dN)[gp].vals, sizeof(double) * info::mesh->dimension * nsize);
			}
			memcpy(w, Mesh::edata[ei->code].weighFactor->data(), sizeof(double) * nsize);
		}

	}

	_invJ.addGeneralInput(_ncoordinates);
	_detJ.addGeneralInput(_ncoordinates);
	_dND.addGeneralInput(_invJ);
	_coordinates.addGeneralInput(_ncoordinates);
	_temperature.addGeneralInput(_ntemperature);
	_isotropicK.addGeneralInput(_conductivity[0]);
	_K.addGeneralInputs(_conductivity[1], _conductivity[2], _conductivity[3], _csCartesian, _csSpherical, _csCylindrical);

	_u.addGeneralInputs(_translationMotion, _dens, _hc);
	_advection.addGeneralInputs(_dND, _u);
	_stiffness.addGeneralInputs(_dND, _isotropicK, _K, _advection);

	eslog::info("\n ============================================================================================= \n");
	eslog::info("  PHYSICS                                                                    HEAT TRANSFER 2D  \n");
	eslog::info(" ============================================================================================= \n");

	if (info::mesh->dimension == 2) {
		setMaterials(info::ecf->heat_transfer_2d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->heat_transfer_2d.material_set);
		validateRegionSettings("THICKNESS", info::ecf->heat_transfer_2d.thickness);
		validateRegionSettings("INITIAL TEMPERATURE", info::ecf->heat_transfer_2d.initial_temperature);
	}
	if (info::mesh->dimension == 3) {
		setMaterials(info::ecf->heat_transfer_3d.material_set);
		validateRegionSettings("MATERIAL", info::ecf->heat_transfer_3d.material_set);
		validateRegionSettings("INITIAL TEMPERATURE", info::ecf->heat_transfer_3d.initial_temperature);
	}

	if (step::loadstep == 0) {
		eslog::info("  MATERIALS                                                                                    \n");
		eslog::info(" --------------------------------------------------------------------------------------------- \n");
		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			eslog::info(" --- %s ---%*s \n", info::mesh->materials[i]->name.c_str(), 84 - info::mesh->materials[i]->name.size(), "");
			const MaterialConfiguration *mat = info::mesh->materials[i];

			switch (mat->coordinate_system.type) {
			case CoordinateSystemConfiguration::TYPE::CARTESIAN:
				eslog::info("    COORDINATE SYSTEM:                                                              CARTESIAN \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, _csCartesian, 0);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "ROTATION.X", mat->coordinate_system.rotation.x, _csCartesian, 0);
					examineMaterialParameter(mat->name, "ROTATION.Y", mat->coordinate_system.rotation.y, _csCartesian, 1);
					examineMaterialParameter(mat->name, "ROTATION.Z", mat->coordinate_system.rotation.z, _csCartesian, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::SPHERICAL:
				if (info::mesh->dimension == 2) {
					eslog::error("HEAT TRANSFER 2D does not support SPHERICAL coordinate system.\n");
				}
				if (info::mesh->dimension == 3) {
					eslog::info("    COORDINATE SYSTEM:                                                              SPHERICAL \n");
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, _csSpherical, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, _csSpherical, 1);
					examineMaterialParameter(mat->name, "CENTER.Z", mat->coordinate_system.center.z, _csSpherical, 2);
				}
				break;
			case CoordinateSystemConfiguration::TYPE::CYLINDRICAL:
				eslog::info("    COORDINATE SYSTEM:                                                            CYLINDRICAL \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, _csCylindrical, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, _csCylindrical, 1);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "CENTER.X", mat->coordinate_system.center.x, _csCylindrical, 0);
					examineMaterialParameter(mat->name, "CENTER.Y", mat->coordinate_system.center.y, _csCylindrical, 1);
				}
				break;
			}
			eslog::info("                                                                                               \n");

			examineMaterialParameter(mat->name, "DENSITY", mat->density, _dens, 0);
			examineMaterialParameter(mat->name, "HEAT CAPACITY", mat->heat_capacity, _hc, 0);
			eslog::info("                                                                                               \n");

			switch (mat->thermal_conductivity.model) {
			case ThermalConductivityConfiguration::MODEL::ISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                              ISOTROPIC \n");
				examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), _conductivity[0], 0);
				break;
			case ThermalConductivityConfiguration::MODEL::DIAGONAL:
				eslog::info("         CONDUCTIVITY:                                                               DIAGONAL \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), _conductivity[1], 0);
					examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), _conductivity[1], 1);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), _conductivity[1], 0);
					examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), _conductivity[1], 1);
					examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), _conductivity[1], 2);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::SYMMETRIC:
				eslog::info("         CONDUCTIVITY:                                                              SYMMETRIC \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), _conductivity[2], 0);
					examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), _conductivity[2], 1);
					examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), _conductivity[2], 2);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), _conductivity[2], 0);
					examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), _conductivity[2], 1);
					examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), _conductivity[2], 2);
					examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), _conductivity[2], 3);
					examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), _conductivity[2], 4);
					examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), _conductivity[2], 5);
				}
				break;
			case ThermalConductivityConfiguration::MODEL::ANISOTROPIC:
				eslog::info("         CONDUCTIVITY:                                                            ANISOTROPIC \n");
				if (info::mesh->dimension == 2) {
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), _conductivity[3], 0);
					examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), _conductivity[3], 1);
					examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), _conductivity[3], 2);
					examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(1, 0), _conductivity[3], 3);
				}
				if (info::mesh->dimension == 3) {
					examineMaterialParameter(mat->name, "KXX", mat->thermal_conductivity.values.get(0, 0), _conductivity[3], 0);
					examineMaterialParameter(mat->name, "KYY", mat->thermal_conductivity.values.get(1, 1), _conductivity[3], 1);
					examineMaterialParameter(mat->name, "KZZ", mat->thermal_conductivity.values.get(2, 2), _conductivity[3], 2);
					examineMaterialParameter(mat->name, "KXY", mat->thermal_conductivity.values.get(0, 1), _conductivity[3], 3);
					examineMaterialParameter(mat->name, "KYZ", mat->thermal_conductivity.values.get(1, 2), _conductivity[3], 4);
					examineMaterialParameter(mat->name, "KXZ", mat->thermal_conductivity.values.get(0, 2), _conductivity[3], 5);
					examineMaterialParameter(mat->name, "KYX", mat->thermal_conductivity.values.get(1, 0), _conductivity[3], 6);
					examineMaterialParameter(mat->name, "KZY", mat->thermal_conductivity.values.get(2, 1), _conductivity[3], 7);
					examineMaterialParameter(mat->name, "KZX", mat->thermal_conductivity.values.get(2, 0), _conductivity[3], 8);
				}

				break;
			}
			eslog::info("                                                                                               \n");
		}

		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		if (info::mesh->dimension == 2) {
			printMaterials(info::ecf->heat_transfer_2d.material_set);
		}
		if (info::mesh->dimension == 3) {
			printMaterials(info::ecf->heat_transfer_3d.material_set);
		}
		eslog::info(" ============================================================================================= \n");

		if (info::mesh->dimension == 2) {
			examineElementParameter("THICKNESS", info::ecf->heat_transfer_2d.thickness, _thickness);
			examineElementParameter("INITIAL TEMPERATURE", info::ecf->heat_transfer_2d.initial_temperature, _ninitTemperature);
			if (configuration.translation_motions.size()) {
				switch (info::ecf->heat_transfer_2d.stabilization) {
				case HeatTransferGlobalSettings::STABILIZATION::CAU: eslog::info("  %s:%77s \n", "STABILIZATION", "CAU"); break;
				case HeatTransferGlobalSettings::STABILIZATION::SUPG: eslog::info("  %s:%77s \n", "STABILIZATION", "SUPG"); break;
				}
				eslog::info("  %s:%85g \n", "SIGMA", info::ecf->heat_transfer_2d.sigma);
			}
		}
		if (info::mesh->dimension == 3) {
			examineElementParameter("INITIAL TEMPERATURE", info::ecf->heat_transfer_3d.initial_temperature, _ninitTemperature);
			if (configuration.translation_motions.size()) {
				switch (info::ecf->heat_transfer_3d.stabilization) {
				case HeatTransferGlobalSettings::STABILIZATION::CAU: eslog::info("  %s:%77s \n", "STABILIZATION", "CAU"); break;
				case HeatTransferGlobalSettings::STABILIZATION::SUPG: eslog::info("  %s:%77s \n", "STABILIZATION", "SUPG"); break;
				}
				eslog::info("  %s:%85g \n", "SIGMA", info::ecf->heat_transfer_3d.sigma);
			}
		}
		_ninitTemperature.swapInput(_coordinates, _ncoordinates);
		eslog::info(" ============================================================================================= \n");
	}

	if (configuration.translation_motions.size()) {
		examineElementParameter("TRANSLATION MOTION.X", configuration.translation_motions, _translationMotion, 0);
		examineElementParameter("TRANSLATION MOTION.Y", configuration.translation_motions, _translationMotion, 1);
		if (info::mesh->dimension == 3) {
			examineElementParameter("TRANSLATION MOTION.Z", configuration.translation_motions, _translationMotion, 2);
		}
	}

	examineBoundaryParameter("TEMPERATURE", configuration.temperature, _btemperature);
	examineBoundaryParameter("HEAT FLOW", configuration.heat_flow, _heatflow);
	examineBoundaryParameter("HEAT FLUX", configuration.heat_flux, _heatflux);

	if (configuration.convection.size()) {
		eslog::info("  CONVECTION                                                                                   \n");
		for (auto it = configuration.convection.begin(); it != configuration.convection.end(); ++it) {
			switch (it->second.type) {
			case ConvectionConfiguration::TYPE::EXTERNAL_FORCED: eslog::info("   %30s:  %57s \n", it->first.c_str(), "EXTERNAL_FORCED"); break;
			case ConvectionConfiguration::TYPE::EXTERNAL_NATURAL: eslog::info("   %30s:  %57s \n", it->first.c_str(), "EXTERNAL_NATURAL"); break;
			case ConvectionConfiguration::TYPE::INTERNAL_FORCED: eslog::info("   %30s:  %57s \n", it->first.c_str(), "INTERNAL_FORCED"); break;
			case ConvectionConfiguration::TYPE::INTERNAL_NATURAL: eslog::info("   %30s:  %57s \n", it->first.c_str(), "INTERNAL_NATURAL"); break;
			case ConvectionConfiguration::TYPE::USER: eslog::info("   %30s:  %57s \n", it->first.c_str(), "USER"); break;
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}
	if (configuration.diffuse_radiation.size()) {
		eslog::info("  DIFFUSE RADIATION                                                                            \n");
		for (auto it = configuration.diffuse_radiation.begin(); it != configuration.diffuse_radiation.end(); ++it) {
			if (it->second.external_temperature.evaluator->parameters.size()) {
				std::string params = Parser::join(", ", it->second.external_temperature.evaluator->parameters);
				eslog::info("   %30s: EX. TEMP. %*s       FNC( %s )\n", it->first.c_str(), 34 - params.size(), "", params.c_str());
			} else {
				eslog::info("   %30s: EX. TEMP. %48g \n", it->first.c_str(), it->second.external_temperature.evaluator->eval(Evaluator::Params()));
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}

	eslog::info(" ============================================================================================= \n");
	eslog::info("  PHYSICS CONFIGURATION VALIDATION                                                       PASS  \n");
	eslog::info(" ============================================================================================= \n");
}

template<typename Ttype>
void HeatTransferKernel::validateRegionSettings(const std::string &name, const std::map<std::string, Ttype> &settings)
{
	for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei) {
		if (ei->region == -1) { // intersected regions
			std::vector<int> set;
			for (auto rindex = ei->regions.begin(); rindex != ei->regions.end(); ++rindex) {
				if (settings.find(info::mesh->elementsRegions[*rindex]->name) != settings.end()) {
					set.push_back(*rindex);
				}
			}
			if (set.size() > 1) {
				std::string regions = info::mesh->elementsRegions[set[0]]->name + ", " + info::mesh->elementsRegions[set[1]]->name;
				eslog::info("  INVALID %s SETTINGS (INTERSECTED REGIONS): %*s\n", name.c_str(), 50 - name.size(), regions.c_str());
				eslog::info(" ============================================================================================= \n");
				eslog::error("ESPRESO error: invalid configuration.\n");
			}
		}
	}
}

void HeatTransferKernel::setMaterials(const std::map<std::string, std::string> &settings)
{
	for (auto ei = info::mesh->elements->eintervals.begin(); ei != info::mesh->elements->eintervals.end(); ++ei) {
		int region = ei->region;
		if (region == -1) { // intersected regions
			for (auto rindex = ei->regions.begin(); rindex != ei->regions.end(); ++rindex) {
				if (settings.find(info::mesh->elementsRegions[*rindex]->name) != settings.end()) {
					region = *rindex;
				}
			}
		}

		if (settings.find(info::mesh->elementsRegions[region]->name) == settings.end()) {
			region = 0;
		}

		auto mat = settings.find(info::mesh->elementsRegions[region]->name);

		if (mat == settings.end()) {
			eslog::error("Invalid material configuration: a region without a material settings found.\n");
		}

		for (size_t i = 0; i < info::mesh->materials.size(); ++i) {
			if (StringCompare::caseInsensitiveEq(info::mesh->materials[i]->name, mat->second)) {
				ei->material = i;
			}
		}
	}
}

void HeatTransferKernel::printMaterials(const std::map<std::string, std::string> &settings)
{
	for (auto reg = info::mesh->elementsRegions.begin() + 1; reg != info::mesh->elementsRegions.end(); ++reg) {
		auto ms = settings.find((*reg)->name);
		if (ms != settings.end()) {
			eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
		} else {
			ms = settings.find("ALL_ELEMENTS");
			if (ms != settings.end()) {
				eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
			} else {
				eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), info::mesh->materials.front()->name.c_str());
			}
		}
	}
}

void HeatTransferKernel::examineMaterialParameter(const std::string &material, const std::string &name, const ECFExpression &settings, ParameterInfo &info, int dimension)
{
	if (settings.evaluator->parameters.size()) {
		std::string params = Parser::join(", ", settings.evaluator->parameters);
		eslog::info("   %18s:  %*s       FNC( %s )\n", name.c_str(), 55 - params.size(), "", params.c_str());
	} else {
		eslog::info("   %18s:  %69g \n", name.c_str(), settings.evaluator->eval(Evaluator::Params()));
	}

	for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
		if (StringCompare::caseInsensitiveEq(info::mesh->materials[info::mesh->elements->eintervals[i].material]->name, material)) {
			info.isset[info.dimensions * i + dimension] = true;
			if (settings.evaluator->parameters.size()) {
				info.isconstant[info.dimensions * i + dimension] = false;
				info.evaluator[info.dimensions * i + dimension] = settings.evaluator;
				insertDependency(info, settings.evaluator->parameters, i, dimension);
			} else {
				info.ivalues[info.dimensions * i + dimension] = settings.evaluator->eval(Evaluator::Params());
			}
		}
	}
}

template<class TSecond>
void HeatTransferKernel::examineElementParameter(const std::string &name, const std::map<std::string, TSecond> &settings, ParameterInfo &info, int dimension, std::function<const Evaluator*(const TSecond &expr)> getevaluator)
{
	if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
		const Evaluator *evaluator = getevaluator(settings.begin()->second);
		if (evaluator->parameters.size()) {
			std::string params = Parser::join(", ", evaluator->parameters);
			eslog::info("  %s:  %*s       FNC( %s )\n", name.c_str(), 65 - params.size(), "", params.c_str());
			info.constness = ParameterInfo::Status::EXPRESSION;
			for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
				info.isset[info.dimensions * i + dimension] = true;
				info.isconstant[info.dimensions * i + dimension] = false;
				info.evaluator[info.dimensions * i + dimension] = evaluator;
				insertDependency(info, evaluator->parameters, i, dimension);
			}
		} else {
			eslog::info("  %s:  %*g \n", name.c_str(), 88 - name.size(), evaluator->eval(Evaluator::Params()));
			info.constness = ParameterInfo::Status::GLOBAL;
			for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
				info.isset[info.dimensions * i + dimension] = true;
				info.isconstant[info.dimensions * i + dimension] = true;
				info.ivalues[info.dimensions * i + dimension] = evaluator->eval(Evaluator::Params());
			}
		}
	} else {
		if (settings.size() == 0) {
			info.constness = ParameterInfo::Status::GLOBAL;
			eslog::info("  %s:  %*g \n", name.c_str(), 88 - name.size(), info.ivalues[dimension]);
		} else {
			eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
			info.constness = ParameterInfo::Status::PER_INTERVAL;
			int rindex = 1, rlast = info::mesh->elementsRegions.size() - 1;
			for (auto reg = info::mesh->elementsRegions.begin() + 1; reg != info::mesh->elementsRegions.end(); ++reg, ++rindex) {
				auto ms = settings.find((*reg)->name);
				if (ms == settings.end()) {
					ms = settings.find("ALL_ELEMENTS");
				}
				if (ms != settings.end()) {
					const Evaluator *evaluator = getevaluator(ms->second);
					if (evaluator->parameters.size()) {
						std::string params = Parser::join(", ", evaluator->parameters);
						eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
					} else {
						eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), evaluator->eval(Evaluator::Params()));
					}
					for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
						if (info::mesh->elements->eintervals[i].region == rindex || (info::mesh->elements->eintervals[i].region == 0 && rindex == rlast)) {
							info.isset[info.dimensions * i + dimension] = true;
							if (evaluator->parameters.size()) {
								info.isconstant[info.dimensions * i + dimension] = false;
								info.evaluator[info.dimensions * i + dimension] = evaluator;
								insertDependency(info, evaluator->parameters, i, dimension);
							} else {
								info.ivalues[info.dimensions * i + dimension] = evaluator->eval(Evaluator::Params());
							}
						}
						if (info::mesh->elements->eintervals[i].region == -1) {
							const std::vector<int> &regions = info::mesh->elements->eintervals[i].regions;
							bool other = false;
							for (auto it = regions.begin(); it != regions.end(); ++it) {
								if (*it != rindex && settings.find((info::mesh->elementsRegions[*it])->name) != settings.end()) {
									other = true;
								}
							}
							if (!other) {
								info.isset[info.dimensions * i + dimension] = true;
								if (evaluator->parameters.size()) {
									info.isconstant[info.dimensions * i + dimension] = false;
									info.evaluator[info.dimensions * i + dimension] = evaluator;
									insertDependency(info, evaluator->parameters, i, dimension);
								} else {
									info.ivalues[info.dimensions * i + dimension] = evaluator->eval(Evaluator::Params());
								}
							}
						}
					}
				} else {
					eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), info.ivalues[info.dimensions + dimension]);
				}
			}
			eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		}
	}
}

void HeatTransferKernel::examineBoundaryParameter(const std::string &name, const std::map<std::string, ECFExpression> &settings, ParameterInfo &info)
{
	if (settings.size()) {
		eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");

		info.ivalues.resize(info::mesh->boundaryRegions.size());
		info.isset.resize(info::mesh->boundaryRegions.size(), false);
		info.isconstant.resize(info::mesh->boundaryRegions.size(), true);
		info.evaluator.resize(info::mesh->boundaryRegions.size(), NULL);

		int rindex = 0;
		for (auto reg = info::mesh->boundaryRegions.begin(); reg != info::mesh->boundaryRegions.end(); ++reg, ++rindex) {
			auto ms = settings.find((*reg)->name);
			if (ms != settings.end()) {
				info.isset[rindex] = true;
				const Evaluator *evaluator = ms->second.evaluator;
				if (evaluator->parameters.size()) {
					std::string params = Parser::join(", ", evaluator->parameters);
					eslog::info("   %30s:  %*s       FNC( %s )\n", (*reg)->name.c_str(), 43 - params.size(), "", params.c_str());
					info.isconstant[rindex] = false;
					info.evaluator[rindex] = evaluator;
					insertDependency(info, evaluator->parameters, rindex);
				} else {
					eslog::info("   %30s:  %57g \n", (*reg)->name.c_str(), evaluator->eval(Evaluator::Params()));
					info.ivalues[rindex] = evaluator->eval(Evaluator::Params());
				}
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
	}
}

void HeatTransferKernel::insertDependency(ParameterInfo &info, const std::vector<std::string> &parameters, esint interval, int dimension)
{
	for (auto it = parameters.begin(); it != parameters.end(); ++it) {
		if (StringCompare::caseInsensitiveEq("INITIAL_TEMPERATURE", *it)) {
			info.addInput(_initTemperature, interval, dimension); break;
		}
		if (StringCompare::caseInsensitiveEq("TEMPERATURE", *it)) {
			info.addInput(_temperature, interval, dimension); break;
		}
		if (StringCompare::caseInsensitiveEq("COORDINATE_X", *it)) {
			info.addInput({ _coordinates, 0, info::mesh->dimension }, interval, dimension); break;
		}
		if (StringCompare::caseInsensitiveEq("COORDINATE_Y", *it)) {
			info.addInput({ _coordinates, 1, info::mesh->dimension }, interval, dimension); break;
		}
		if (StringCompare::caseInsensitiveEq("COORDINATE_Z", *it)) {
			info.addInput({ _coordinates, 2, info::mesh->dimension }, interval, dimension); break;
		}
		if (StringCompare::caseInsensitiveEq("TIME", *it)) {
			info.addInput(_time, interval, dimension); break;
		}
		eslog::error("ESPRESO internal error: implement dependency on parameter: '%s'\n", it->c_str());
	}
}

HeatTransferKernel::~HeatTransferKernel()
{

}
