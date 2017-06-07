
#include "parser.h"
#include "../../mesh/structures/elementtypes.h"
#include "../../mesh/structures/mesh.h"
#include "../../mesh/structures/coordinates.h"
#include "../../mesh/structures/material.h"
#include "../../mesh/structures/region.h"
#include "../../mesh/settings/evaluator.h"
#include "../../mesh/elements/element.h"

#include "../../basis/logging/logging.h"

using namespace espreso::input;

WorkbenchParser::WorkbenchParser(Mesh &mesh): bodyCounter(0), _mesh(mesh)
{
	_commands["/wb"] = WorkbenchCommands::WB;
	_commands["nblock"] = WorkbenchCommands::NBLOCK;
	_commands["eblock"] = WorkbenchCommands::EBLOCK;
	_commands["cmblock"] = WorkbenchCommands::CMBLOCK;
	_commands["MP"] = WorkbenchCommands::MP;
	_commands["MPTEMP"] = WorkbenchCommands::MPTEMP;

	_commands["et"] = WorkbenchCommands::ET;
	_commands["cmsel"] = WorkbenchCommands::CMSEL;
	_commands["*DIM"] = WorkbenchCommands::DIM;
	_commands["d"] = WorkbenchCommands::DIRICHLET;
	_commands["tunif"] = WorkbenchCommands::INITIAL_TEMPERATURE;
	_commands["f"] = WorkbenchCommands::FORCE;
	_commands["sf"] = WorkbenchCommands::SURFACE_EFFECT;
	_commands["acel"] = WorkbenchCommands::ACCELERATION;
	_commands["obstacle"] = WorkbenchCommands::OBSTACLE;
	_commands["nsel"] = WorkbenchCommands::NSEL;
	_commands["esel"] = WorkbenchCommands::ESEL;
	_commands["_loadvari"] = WorkbenchCommands::LOADVAR;
}

WorkbenchParser::~WorkbenchParser()
{
	for (size_t t = 0; t < _tables.size(); t++) {
		delete _tables[t];
	}
}

void WorkbenchParser::open(std::string path)
{
	_file.open(path.c_str());
	if (!_file.is_open()) {
		ESINFO(GLOBAL_ERROR) << "Cannot load mesh from file: " << path;
	}
}

std::vector<std::string> WorkbenchParser::divide(std::string &line, std::string delim)
{
	std::vector<std::string> parts;
	size_t start, end = -1;
	do {
		start = end + 1;
		end = line.find_first_of(delim, start);
		parts.push_back(Parser::strip(line.substr(start, end - start)));
		if (parts.back().size() && parts.back()[0] == '!') {
			parts.pop_back();
			break;
		}
	} while(end != std::string::npos);
	if (!parts.back().size()) {
		parts.pop_back();
	}
	return parts;
}

std::vector<int> WorkbenchParser::parseBlockHeader(std::string &line)
{
	std::vector<int> sizes;
	size_t start, end = 0;
	std::string substr;
	do {
		start = end + 1;
		end = line.find_first_of(",)", start);
		substr = line.substr(start, end - start);
		std::vector<std::string> parts = divide(substr, "ie");
		sizes.insert(sizes.end(), std::stoi(parts[0]), std::stoi(parts[1]));
	} while(end < line.size() && line[end] != ')');

	return sizes;
}

bool WorkbenchParser::trash(const std::string &line)
{
	if (line.size() && line[0] == '!') { // comment
		return false;
	}
	return false;
}

WorkbenchCommands WorkbenchParser::process()
{
	while (_file.good()) {
		do {
			getline(_file, _line);
		} while (trash(_line));

		std::string command = _line.substr(0, _line.find_first_of(','));
		auto it = _commands.find(command);
		if (it != _commands.end()) {
			switch (it->second) {
			case WorkbenchCommands::ET: et(); break;
			case WorkbenchCommands::CMSEL: cmsel(); break;
			case WorkbenchCommands::NSEL: nsel(); break;
			case WorkbenchCommands::ESEL: esel(); break;
			case WorkbenchCommands::DIM: dim(); break;

			default: {
				return it->second;
			}
			}
		}

	}

	return WorkbenchCommands::END;
}

bool WorkbenchParser::workbench(const std::string type, const std::string status)
{
	std::vector<std::string> params = divide(_line);
	return params[1].compare(0, type.size(), type) == 0 && params[2].compare(0, status.size(), status) == 0;
}

void WorkbenchParser::nblock(Coordinates &coordinates)
{
	std::vector<std::string> params = divide(_line);
//	eslocal NUMFILED, Solkey, NDMAX, NDSEL; // TODO: use all parameters
//	switch (params.size()) {
//	case 5:
//		NDSEL = std::stol(params[3]);
//	case 4:
//		NDMAX = std::stol(params[3]);
//	case 3:
//		Solkey = std::stol(params[2]);
//		ESINFO(GLOBAL_ERROR) << "The input point format is not implemented in ESPRESO";
//	case 2:
//		NUMFILED = std::stol(params[1]);
//	}

	getline(_file, _line);
	std::vector<int> sizes = parseBlockHeader(_line);
	size_t offset;

	Point point;
	eslocal id;

	ESINFO(DETAILS) << "WB: CREATE NBLOCK";

	while (true) {
		getline(_file, _line);
		id = std::stol(_line.substr(0, sizes[0]));
		if (id == -1) {
			return; // TODO: loading more NBLOCKs
		}
		offset = sizes[0];

		point.x = std::stod(_line.substr(offset, sizes[1]));
		offset += sizes[1];
		point.y = std::stod(_line.substr(offset, sizes[2]));
		offset += sizes[2];
		point.z = std::stod(_line.substr(offset, sizes[3]));
		offset += sizes[3];

		coordinates.add(point, id - 1, id - 1);
	}
}

void WorkbenchParser::eblock(std::vector<Element*> &elements, std::vector<Region*> &regions, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	regions.push_back(new Region(ElementType::ELEMENTS));

	int MATERIAL = 0, ETYPE = 1, CONSTANT = 2, COORDINATES, NODES, PARAM_SIZE = 6, NODE_SIZE;
	bool SOLID = true;

	std::vector<std::string> params = divide(_line);
	eslocal NUMNODES = 0, NDSEL = 0; //, NDMAX;
	switch (params.size()) {
	case 5:
		NDSEL = params[4].size() ? std::stol(params[4]) : 0;
	case 4:
		//NDMAX = params[3].size() ? std::stol(params[3]) : 0;
	case 3:
		if (!StringCompare::caseInsensitiveEq(params[2], "SOLID")) {
			SOLID = false;
			MATERIAL = 3;
			COORDINATES = 4;
			PARAM_SIZE = 5;
		} else { // SOLID element
			SOLID = true;
			MATERIAL = 0;
			COORDINATES = 4;
			NODES = 8;
			PARAM_SIZE = 11;
		}
		NUMNODES = std::stoi(params[1]);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Unknown eblock format";
	}

	getline(_file, _line);
	std::vector<int> sizes = parseBlockHeader(_line);

	if (NDSEL) {
		elements.reserve(elements.size() + NDSEL);
	}

	ESINFO(DETAILS) << "WB: CREATE EBLOCK";

	std::vector<eslocal> values(38), eParams(Element::PARAMS_SIZE);
	for (eslocal i = 0; i < (NDSEL ? NDSEL : i + 1); i++) {
		getline(_file, _line);
		size_t start = 0;
		NODE_SIZE = 0;
		for (size_t i = 0; i < sizes.size() && _line.size() > start + sizes[i]; i++) {
			values[i] = std::stol(_line.substr(start, sizes[i]));
			if (i >= (size_t)PARAM_SIZE) {
				values[i]--;
				NODE_SIZE++;
			}
			if (NDSEL == 0 && values[i] == -2) { // -2 == -1 - 1
				break;
			}
			start += sizes[i];
		}
		if (SOLID && values[NODES] > NUMNODES - PARAM_SIZE) {
			start = 0;
			getline(_file, _line);
			for (size_t i = 0; i < values[NODES] - sizes.size() + PARAM_SIZE; i++) {
				values[sizes.size() + i] = std::stol(_line.substr(start, sizes[i])) - 1;
				start += sizes[i];
			}
		}

		if (SOLID) {
			NODE_SIZE = values[NODES];
		}

		eParams[Element::MATERIAL] = values[MATERIAL] - 1;
		eParams[Element::CONSTANT] = values[CONSTANT];
		eParams[Element::COORDINATES] = values[COORDINATES] - 1;
		eParams[Element::BODY] = 0; // FIX ME! bodyCounter;
		Element *e = AnsysUtils::createElement(values.data() + PARAM_SIZE, NODE_SIZE, eParams.data(), eType[values[ETYPE] - 1]);
		switch (eType[values[ETYPE] - 1]) {
		case 77:
		case 87:
		case 90:
		case 185:
		case 186:
		case 187:
			elements.push_back(e);
			break;
		case 152:
		case 154:
			faces.push_back(e);
			break;
		case 151:
			edges.push_back(e);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Unknown element type " << eType[values[ETYPE] - 1];
		}
		regions.back()->elements().push_back(e);
	}
	switch (regions.back()->elements().back()->type()) {
	case Element::Type::VOLUME:
		regions.back()->eType = ElementType::ELEMENTS;
		break;
	case Element::Type::PLANE:
		regions.back()->eType = ElementType::FACES;
		break;
	case Element::Type::LINE:
		regions.back()->eType = ElementType::EDGES;
		break;
	default:
		break;
	}

	regions.back()->name = std::to_string(values[ETYPE]);
	bodyCounter++;
}

void WorkbenchParser::mp(std::vector<Material*> &materials, Evaluator *evaluator)
{
	std::vector<std::string> params = divide(_line);

	size_t mNumber = std::stoi(params[2]);
	for (size_t i = materials.size(); i < mNumber; i++) {
		materials.push_back(new Material(_mesh.coordinates()));
		materials.back()->setModel(PHYSICS::ADVECTION_DIFFUSION_2D, MATERIAL_MODEL::ISOTROPIC);
		materials.back()->setModel(PHYSICS::ADVECTION_DIFFUSION_3D, MATERIAL_MODEL::ISOTROPIC);
		materials.back()->setModel(PHYSICS::LINEAR_ELASTICITY_2D, MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC);
		materials.back()->setModel(PHYSICS::LINEAR_ELASTICITY_3D, MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC);
	}
	mNumber--;

	bool match = false;
	auto set = [&] (MATERIAL_PARAMETER parameter, const std::string &name) {
		if (StringCompare::caseInsensitiveEq(params[1], name)) {
			if (evaluator == NULL) {
				materials[mNumber]->set(parameter, params[3]);
			} else {
				materials[mNumber]->set(parameter, evaluator);
			}
			match = true;
			return true;
		}
		return false;
	};

	auto setWithModel = [&] (MATERIAL_PARAMETER parameter, const std::string &name, std::function<void(void)> setModel) {
		if (set(parameter, name)) {
			setModel();
		}
	};

	auto setModelAD = [&] (MATERIAL_MODEL model) {
		materials[mNumber]->setModel(PHYSICS::ADVECTION_DIFFUSION_2D, model);
		materials[mNumber]->setModel(PHYSICS::ADVECTION_DIFFUSION_3D, model);
	};

	auto setModelLE = [&] (MATERIAL_MODEL model) {
		materials[mNumber]->setModel(PHYSICS::LINEAR_ELASTICITY_2D, model);
		materials[mNumber]->setModel(PHYSICS::LINEAR_ELASTICITY_3D, model);
	};

	auto skip = [&] (const std::string &name) {
		if (StringCompare::caseInsensitiveEq(params[1], name)) {
			match = true;
		}
	};

	set(MATERIAL_PARAMETER::DENSITY                , "DENS");
	set(MATERIAL_PARAMETER::HEAT_CAPACITY          , "C");

	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XX, "KXX", [&] () { setModelAD(MATERIAL_MODEL::ISOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YY, "KYY", [&] () { setModelAD(MATERIAL_MODEL::DIAGONAL); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZZ, "KZZ", [&] () { setModelAD(MATERIAL_MODEL::DIAGONAL); });

	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XY, "KXY", [&] () { setModelAD(MATERIAL_MODEL::SYMMETRIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_XZ, "KXZ", [&] () { setModelAD(MATERIAL_MODEL::SYMMETRIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YZ, "KYZ", [&] () { setModelAD(MATERIAL_MODEL::SYMMETRIC); });

	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_YX, "KYX", [&] () { setModelAD(MATERIAL_MODEL::ANISOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZX, "KZX", [&] () { setModelAD(MATERIAL_MODEL::ANISOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_CONDUCTIVITY_ZY, "KZY", [&] () { setModelAD(MATERIAL_MODEL::ANISOTROPIC); });


	setWithModel(MATERIAL_PARAMETER::YOUNG_MODULUS_X        , "EX"  , [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::POISSON_RATIO_XY       , "NUXY", [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_EXPANSION_X    , "ALPX", [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ISOTROPIC); });

	setWithModel(MATERIAL_PARAMETER::YOUNG_MODULUS_Y        , "EY"  , [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::YOUNG_MODULUS_Z        , "EZ"  , [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::POISSON_RATIO_XZ       , "NUXZ", [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::POISSON_RATIO_YZ       , "NUYZ", [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_EXPANSION_Y    , "ALPY", [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC); });
	setWithModel(MATERIAL_PARAMETER::THERMAL_EXPANSION_Z    , "ALPZ", [&] () { setModelLE(MATERIAL_MODEL::LINEAR_ELASTIC_ORTHOTROPIC); });

	skip("RSVX");
	skip("RSVY");
	skip("RSVZ");
	skip("MURX");
	skip("MURY");
	skip("MURZ");

	if (!match) {
		ESINFO(GLOBAL_ERROR) << "Unknown material property '" << params[1] << "'";
	}
}

void WorkbenchParser::mptemp(std::vector<Material*> &materials)
{
	std::vector<std::pair<double, double> > table;
	while (true) {
		std::vector<std::string> params = divide(_line);
		if (!params[0].compare("MPDATA")) {
			size_t loaded = 0;
			while (loaded < table.size()) {
				if (loaded) {
					getline(_file, _line);
					params = divide(_line);
				}
				for (size_t i = 0; i + 4 < params.size(); i++) {
					table[loaded++].second = std::stod(params[i + 4]);
				}
			}
			mp(materials, new TableInterpolationEvaluator("TEMP" + params[2], table));
		} else {
			if (!params[1].size()) {
				break;
			}
			int index = std::stoi(params[1]);
			table.resize(index);
			table[index - 1].first = std::stod(params[2]);
		}
		getline(_file, _line);
	}
}

static size_t getRegionIndex(std::vector<espreso::Region*> &regions, const std::string &name)
{
	if (name.size()) {
		for (size_t i = 0; i < regions.size(); i++) {
			if (espreso::StringCompare::caseInsensitiveEq(regions[i]->name, name)) {
				return i;
			}
		}
	}
	ESINFO(espreso::GLOBAL_ERROR) << "Unknown region '" << name << "'";
	return -1;
}

void WorkbenchParser::cmblock(std::vector<Element*> &elements, std::vector<Region*> &regions, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);

	params[1] = Parser::strip(params[1]);
	regions.push_back(new Region(ElementType::NODES));
	regions.back()->name = params[1];
	std::vector<Element*> &region = regions.back()->elements();
	if (params[2].compare(0, 4, "NODE") == 0) {
		ESINFO(DETAILS) << "WB: CREATE BLOCK OF NODES: " << params[1];
		eslocal size = std::stol(params[3]);

		getline(_file, _line);
		std::vector<int> sizes = parseBlockHeader(_line);

		getline(_file, _line);
		int start = 0, n = 0, number = 0;
		while (number++ < size) {
			region.push_back(nodes[std::stol(_line.substr(start, sizes[n])) - 1]);
			start += sizes[n++];
			if (n % sizes.size() == 0) {
				getline(_file, _line);
				start = 0;
				n = 0;
			}
		}
	} else {
		ESINFO(GLOBAL_ERROR) << "Not implemented loading of cmblock of elements";
	}
}

static void pushTableEvaluator(std::vector<espreso::Evaluator*> &evaluators, const std::string &name, const std::vector<espreso::TableEvaluator*> &tables)
{
	for (size_t t = 0; t < tables.size(); t++) {
		if (!name.compare(1, tables[t]->name().size(), tables[t]->name())) {
			evaluators.push_back(tables[t]->copy());
			return;
		}
	}
	ESINFO(espreso::GLOBAL_ERROR) << "Unknown table '" << name << "'";
}

static void pushEvaluator(std::vector<espreso::Evaluator*> &evaluators, const std::string &value, const espreso::Coordinates &coordinates) {
	if (value.find("x") == std::string::npos && value.find("y") == std::string::npos && value.find("z") == std::string::npos) {
		espreso::Expression expr(value, {});
		evaluators.push_back(new espreso::ConstEvaluator(expr.evaluate({})));
	} else {
		evaluators.push_back(new espreso::CoordinatesEvaluator(value, coordinates));
	}
};

bool WorkbenchParser::setProperty(const std::string &parameter, const std::string &value, const std::string &name, espreso::Property property, size_t loadStep, espreso::Region *region, std::vector<Evaluator*> &evaluators) {
	if (!parameter.compare(0, name.size(), name)) {
		if (!value.compare(0, 1, "%") && !value.compare(value.size() - 1, value.size(), "%")) {
			pushTableEvaluator(evaluators, value, _tables);
			ESINFO(espreso::DETAILS) << "WB: SET " << property << " to region " << region->name << " to table " << evaluators.back()->name();
		} else {
			pushEvaluator(evaluators, value, _mesh.coordinates());
			ESINFO(espreso::DETAILS) << "WB: SET " << property << " to region " << region->name << " to " << value;
		}
		evaluators.back()->property() = property;
		region->settings.resize(loadStep + 1);
		region->settings[loadStep][property].push_back(evaluators.back());
		return true;
	}
	return false;
};

void WorkbenchParser::dirichlet(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);
	bool set = false;

	size_t loadStep = 0;
	if (!params[1].compare(0, 3, "all")) {
		Region *region = regions[getRegionIndex(regions, _selectedRegion.size() ? _selectedRegion : "ALL_NODES")];
		region->settings.resize(loadStep + 1);
		if (!params[2].compare(0, 3, "all")) {
			evaluators.push_back(new ConstEvaluator(0, Property::DISPLACEMENT_X));
			region->settings[loadStep][Property::DISPLACEMENT_X].push_back(evaluators.back());
			region->settings[loadStep][Property::DISPLACEMENT_Y].push_back(evaluators.back());
			region->settings[loadStep][Property::DISPLACEMENT_Z].push_back(evaluators.back());
			return;
		} else {
			set = set || setProperty(params[2], params[3], "ux", Property::DISPLACEMENT_X, loadStep, region, evaluators);
			set = set || setProperty(params[2], params[3], "uy", Property::DISPLACEMENT_Y, loadStep, region, evaluators);
			set = set || setProperty(params[2], params[3], "uz", Property::DISPLACEMENT_Z, loadStep, region, evaluators);
		}
	}

	Region *region = regions[getRegionIndex(regions, params[1])];

	set = set || setProperty(params[2], params[3], "ux"  , Property::DISPLACEMENT_X, loadStep, region, evaluators);
	set = set || setProperty(params[2], params[3], "uy"  , Property::DISPLACEMENT_X, loadStep, region, evaluators);
	set = set || setProperty(params[2], params[3], "uz"  , Property::DISPLACEMENT_X, loadStep, region, evaluators);
	set = set || setProperty(params[2], params[3], "temp", Property::TEMPERATURE   , loadStep, region, evaluators);

	if (!set) {
		ESINFO(GLOBAL_ERROR) << "Not implemented d option";
	}
}

void WorkbenchParser::force(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);
	bool set = false;

	if (!params[1].compare(0, 3, "all")) {
		if (!params[2].compare(0, 3, "all")) {
			ESINFO(GLOBAL_ERROR) << "Broken WORKBENCH INPUT";
		} else {
			ESINFO(GLOBAL_ERROR) << "Not implemented f option '" << params[2] << "'";
		}
	}

	size_t loadStep = 0;
	Region *region = regions[getRegionIndex(regions, params[1])];

	set = set || setProperty(params[2], params[3], "fx", Property::FORCE_X, loadStep, region, evaluators);
	set = set || setProperty(params[2], params[3], "fy", Property::FORCE_Y, loadStep, region, evaluators);
	set = set || setProperty(params[2], params[3], "fz", Property::FORCE_Z, loadStep, region, evaluators);

	if (!set) {
		ESINFO(GLOBAL_ERROR) << "Not implemented d option '" << params[2] << "'";
	}
}

void WorkbenchParser::sf(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	std::vector<std::string> params = divide(_line);
	bool set = false;

	size_t loadStep = 0;
	Region *region = regions[getRegionIndex(regions, _selectedRegion)];

	if (!params[1].compare(0, 3, "all")) {
		set = set || setProperty(params[2], params[3], "pres", Property::PRESSURE, loadStep, region, evaluators);
		set = set || setProperty(params[2], params[3], "hflux", Property::HEAT_FLUX, loadStep, region, evaluators);
		if (!params[2].compare(0, 4, "conv")) {
			set = true;
			set = set && setProperty(params[2], params[3], "conv", Property::HEAT_TRANSFER_COEFFICIENT, loadStep, region, evaluators);
			set = set && setProperty(params[2], params[4], "conv", Property::EXTERNAL_TEMPERATURE, loadStep, region, evaluators);
		}
	}

	if (!set) {
		ESINFO(GLOBAL_ERROR) << "Not implemented sf option '" << params[2] << "'";
	}
}

void WorkbenchParser::acceleration(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions)
{
	std::vector<std::string> params = divide(_line);
	bool set = true;

	size_t loadStep = 0;
	Region *region = regions[getRegionIndex(regions, _selectedRegion.size() ? _selectedRegion : "ALL_ELEMENTS")];

	set = set && setProperty(params[0], params[1], "acel", Property::ACCELERATION_X, loadStep, region, evaluators);
	set = set && setProperty(params[0], params[2], "acel", Property::ACCELERATION_Y, loadStep, region, evaluators);
	set = set && setProperty(params[0], params[3], "acel", Property::ACCELERATION_Z, loadStep, region, evaluators);

	if (!set) {
		ESINFO(GLOBAL_ERROR) << "Unexpected acel option";
	}
}

void WorkbenchParser::initial_temperature(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions)
{
	std::vector<std::string> params = divide(_line);
	bool set = true;

	size_t loadStep = 0;
	Region *region = regions[getRegionIndex(regions, _selectedRegion.size() ? _selectedRegion : "ALL_ELEMENTS")];

	set = set && setProperty(params[0], params[1], "tunif", Property::INITIAL_TEMPERATURE, loadStep, region, evaluators);

	if (!set) {
		ESINFO(GLOBAL_ERROR) << "Unexpected tunif option";
	}
}

void WorkbenchParser::obstacle(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);

	if (!_selectedRegion.size()) {
		ESINFO(GLOBAL_ERROR) << "You have to set region for obstacle";
	}

	size_t loadStep = 0;
	Region *region = regions[getRegionIndex(regions, _selectedRegion.size() ? _selectedRegion : "ALL_NODES")];
	region->settings.resize(loadStep + 1);
	evaluators.push_back(new CoordinatesEvaluator(params[1], _mesh.coordinates(), Property::OBSTACLE));
	evaluators.push_back(new CoordinatesEvaluator(params[2], _mesh.coordinates(), Property::NORMAL_DIRECTION));

	ESINFO(DETAILS) << "WB: SET obstacle to region " << _selectedRegion << " to " << params[1] << ", direction: " << params[2];

	region->settings[loadStep][Property::OBSTACLE].push_back(evaluators[evaluators.size() - 2]);
	region->settings[loadStep][Property::NORMAL_DIRECTION].push_back(evaluators[evaluators.size() - 1]);
}

void WorkbenchParser::loadvar()
{

}

void WorkbenchParser::dim()
{
	std::vector<std::string> params = divide(_line);

	if (!params[2].compare(0, 5, "TABLE")) {
		ESINFO(DETAILS) << "WB: CREATE TABLE: " << params[1];
		std::vector<std::vector<std::vector<double> > > table;
		std::vector<TableEvaluator::TableProperty> properties;
		std::vector<std::vector<double> > axis;
		if (3 < params.size()) {
			table.resize(std::stoi(params[3]));
		}
		if (4 < params.size()) {
			for (size_t i = 0; i < table.size(); i++) {
				table[i].resize(std::stoi(params[4]));
			}
		}
		if (5 < params.size()) {
			for (size_t i = 0; i < table.size(); i++) {
				for (size_t j = 0; j < table[i].size(); j++) {
					table[i][j].resize(std::stoi(params[5]));
				}
			}
		}
		for (size_t p = 6; p < params.size(); p++) {
			if (!params[p].compare(0, 4, "TIME")) {
				properties.push_back(TableEvaluator::TableProperty::TIME);
			}
			if (!params[p].compare(0, 4, "TEMP")) {
				properties.push_back(TableEvaluator::TableProperty::TEMPERATURE);
			}
			if (!params[p].compare(0, 8, "PRESSURE")) {
				properties.push_back(TableEvaluator::TableProperty::PRESSURE);
			}
			if (!params[p].compare(0, 8, "VELOCITY")) {
				properties.push_back(TableEvaluator::TableProperty::VELOCITY);
			}
		}
		size_t setted = 0;
		while (setted < table.size() * table[0].size() * table[0][0].size()) {
			getline(_file, _line);
			if (!_line.compare(0, params[1].size(), params[1])) {
				std::string indices = _line.substr(_line.find_first_of("(") + 1, _line.find_first_of(")") - _line.find_first_of("(") - 1);
				std::string value = _line.substr(_line.find_first_of("=") + 1);
				std::vector<std::string> parsed_indices = divide(indices);
				std::vector<size_t> index;
				for (size_t i = 0; i < parsed_indices.size(); i++) {
					index.push_back(std::stoi(parsed_indices[i]));
				}
				if (index[1] < 1) {
					axis.resize(1, std::vector<double>(table.size()));
					axis[0][index[0] - 1] = std::stod(value);
				}
				if (index[2] < 1) {
					axis.resize(2, std::vector<double>(table[0].size()));
					axis[1][index[1] - 1] = std::stod(value);
				}

				if (index[1] > 0 && index[2] > 0) {
					setted++;
					table[index[0] - 1][index[1] - 1][index[2] - 1] = std::stod(value);
				}
			}
		}

		_tables.push_back(new TableEvaluator(params[1], table, properties, axis));
		return;
	}

	if (!params[2].compare(0, 5, "string")) {
		ESINFO(ALWAYS) << Info::TextColor::YELLOW << "Skip *DIM option 'string'";
		return;
	}
}

void WorkbenchParser::et()
{
	std::vector<std::string> params = divide(_line);
	int n = std::stoi(params[1]);
	eType.resize(n);
	eType[n - 1] = std::stoi(params[2]);
	ESINFO(DETAILS) << "WB: SET ELEMENT TYPE: " << eType[n - 1];
}

void WorkbenchParser::cmsel()
{
	std::vector<std::string> params = divide(_line);
	if (params[1].compare(0, 1, "s")) {
		ESINFO(GLOBAL_ERROR) << Info::TextColor::YELLOW << "unknown cmsel parameters '" << params[1] << "'";
	}
	_selectedRegion = Parser::strip(Parser::split(params[2], "!")[0]);
	ESINFO(DETAILS) << "WB: SELECT REGION: " << _selectedRegion;
}

void WorkbenchParser::nsel()
{
	std::vector<std::string> params = divide(_line);

	if (params[1].compare(0, 3, "all") == 0) {
		_selectedRegion.clear();
		ESINFO(DETAILS) << "WB: CLEAR SELECTION";
	} else {
		ESINFO(GLOBAL_ERROR) << "Not implemented nsel option '" << params[1] << "'";
	}
}

void WorkbenchParser::esel()
{
	std::vector<std::string> params = divide(_line);

	if (params[1].compare(0, 3, "all") == 0) {
		_selectedRegion.clear();
		ESINFO(DETAILS) << "WB: CLEAR SELECTION";
		return;
	}
	if (params[1].compare(0, 4, "none") == 0) {
		_selectedRegion.clear();
		ESINFO(DETAILS) << "WB: CLEAR SELECTION";
		return;
	}
	if (params[1].compare(0, 1, "a") == 0) {
		ESINFO(OVERVIEW) << Info::TextColor::YELLOW << "Skip selection of material: " << _line;
		return;
	}
	if (params[1].compare(0, 1, "s") == 0 && params[2].compare(0, 3, "mat") == 0) {
		ESINFO(OVERVIEW) << Info::TextColor::YELLOW << "Skip selection of material: " << _line;
		return;
	}
	if (params[1].compare(0, 1, "s") == 0 && params[2].compare(0, 4, "type") == 0) {
		_selectedRegion = params[4];
		ESINFO(DETAILS) << "WB: SELECT ELEMENTS REGION " << params[4];
		return;
	}

	ESINFO(GLOBAL_ERROR) << "Not implemented esel option '" << _line << "'";
}




