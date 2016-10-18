
#include "parser.h"

using namespace espreso::input;

WorkbenchParser::WorkbenchParser(Mesh &mesh): bodyCounter(0), _mesh(mesh)
{
	_commands["/wb"] = WorkbenchCommands::WB;
	_commands["nblock"] = WorkbenchCommands::NBLOCK;
	_commands["eblock"] = WorkbenchCommands::EBLOCK;
	_commands["cmblock"] = WorkbenchCommands::CMBLOCK;
	_commands["MP"] = WorkbenchCommands::MP;

	_commands["et"] = WorkbenchCommands::ET;
	_commands["cmsel"] = WorkbenchCommands::CMSEL;
	_commands["*DIM"] = WorkbenchCommands::DIM;
	_commands["d"] = WorkbenchCommands::DISPLACEMENT;
	_commands["f"] = WorkbenchCommands::FORCE;
	_commands["obstacle"] = WorkbenchCommands::OBSTACLE;
	_commands["nsel"] = WorkbenchCommands::NSEL;
	_commands["_loadvari"] = WorkbenchCommands::LOADVAR;
}

std::vector<std::string> WorkbenchParser::divide(std::string &line, std::string delim)
{
	std::vector<std::string> parts;
	size_t start, end = -1;
	do {
		start = end + 1;
		end = line.find_first_of(delim, start);
		parts.push_back(Parser::strip(line.substr(start, end - start)));
	} while(end != std::string::npos);

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
	eslocal NUMFILED, Solkey, NDMAX, NDSEL; // TODO: use all parameters
	switch (params.size()) {
	case 5:
		NDSEL = std::stol(params[3]);
	case 4:
		NDMAX = std::stol(params[3]);
	case 3:
		Solkey = std::stol(params[2]);
		ESINFO(GLOBAL_ERROR) << "The input point format is not implemented in ESPRESO";
	case 2:
		NUMFILED = std::stol(params[1]);
	}

	getline(_file, _line);
	std::vector<int> sizes = parseBlockHeader(_line);
	size_t offset;

	Point point;
	size_t id;

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

void WorkbenchParser::eblock(std::vector<Element*> &elements)
{
	int MATERIAL, ETYPE, CONSTANT, COORDINATES, NODES, PARAM_SIZE, NODE_SIZE;
	bool SOLID;

	std::vector<std::string> params = divide(_line);
	eslocal NUMNODES, NDMAX, NDSEL;
	switch (params.size()) {
	case 5:
		NDSEL = params[4].size() ? std::stol(params[4]) : 0;
	case 4:
		NDMAX = params[3].size() ? std::stol(params[3]) : 0;
	case 3:
		if (!StringCompare::caseInsensitiveEq(params[2], "SOLID")) {
			SOLID = false;
			ETYPE = 1;
			CONSTANT = 2;
			MATERIAL = 3;
			COORDINATES = 4;
			PARAM_SIZE = 5;
		} else { // SOLID element
			SOLID = true;
			MATERIAL = 0;
			ETYPE = 1;
			CONSTANT = 2;
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
		int start = 0;
		NODE_SIZE = 0;
		for (size_t i = 0; i < sizes.size() && _line.size() > start + sizes[i]; i++) {
			values[i] = std::stol(_line.substr(start, sizes[i]));
			if (i >= PARAM_SIZE) {
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
		eParams[Element::BODY] = bodyCounter;
		elements.push_back(AnsysUtils::createElement(values.data() + PARAM_SIZE, NODE_SIZE, eParams.data(), eType[values[ETYPE] - 1]));
	}
	bodyCounter++;
}

void WorkbenchParser::mp(std::vector<Material> &materials)
{
	std::vector<std::string> params = divide(_line);

	int mNumber = std::stoi(params[2]);
	materials.resize(mNumber--, Material(_mesh.coordinates()));

	if (!materials[mNumber].setParameter(params[1], params[3])) {
		ESINFO(GLOBAL_ERROR) << "Unknown material property '" << params[1] << "'";
	}
}

void WorkbenchParser::eblock(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	ESINFO(GLOBAL_ERROR) << "Implement eblock settings";
	//SurfaceCondition *ec = new SurfaceCondition();
	//eblock(ec->faces());

	//conditions.push_back(ec);
//	std::stringstream ss;
//	ss << eType.back();
//	selections.push_back({ss.str(), ConditionElements::ELEMENTS});
}

void WorkbenchParser::cmblock(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);

	params[1] = Parser::strip(params[1]);
	std::vector<Element*> &region = _regions[params[1]];
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

static void pushTableEvaluator(std::vector<espreso::Evaluator*> &evaluators, const std::vector<espreso::TableEvaluator*> &tables, const std::string &name)
{
	for (size_t t = 0; t < tables.size(); t++) {
		if (!name.compare(1, tables[t]->name().size(), tables[t]->name())) {
			evaluators.push_back(tables[t]->copy());
			return;
		}
	}
	ESINFO(espreso::GLOBAL_ERROR) << "Unknown table '" << name << "'";
}

void WorkbenchParser::displacement(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);

	auto pushEvaluator = [&] (const std::string &value, Property property) {
		if (value.find("x") == std::string::npos && value.find("y") == std::string::npos && value.find("z") == std::string::npos) {
			espreso::Expression expr(value, {});
			evaluators.push_back(new espreso::ConstEvaluator(expr.evaluate({}), property));
		} else {
			evaluators.push_back(new espreso::CoordinatesEvaluator(value, _mesh.coordinates(), property));
		}
	};

	if (!params[1].compare(0, 3, "all")) {
		std::vector<Element*> &region = _selectedRegion.size() ? _regions[_selectedRegion] : nodes;
		if (!params[2].compare(0, 3, "all")) {
			evaluators.push_back(new ConstEvaluator(0, Property::DISPLACEMENT_X)); // TODO: make displacement property
			for (size_t e = 0; e < region.size(); e++) {
				region[e]->addSettings(Property::DISPLACEMENT_X, evaluators.back());
				region[e]->addSettings(Property::DISPLACEMENT_Y, evaluators.back());
				region[e]->addSettings(Property::DISPLACEMENT_Z, evaluators.back());
			}
			return;
		} else {
			if (params[3].compare(0, 1, "%") || params[3].compare(params[3].size() - 1, params[3].size(), "%")) {
				if (!params[2].compare(0, 2, "ux")) {
					ESINFO(DETAILS) << "WB: SET ux to " << (_selectedRegion.size() ? "region " + _selectedRegion : "all nodes") << " to " << params[3];
					pushEvaluator(params[3], Property::DISPLACEMENT_X);
					for (size_t e = 0; e < region.size(); e++) {
						region[e]->addSettings(Property::DISPLACEMENT_X, evaluators.back());
					}
					return;
				}
				if (!params[2].compare(0, 2, "uy")) {
					ESINFO(DETAILS) << "WB: SET uy to " << (_selectedRegion.size() ? "region " + _selectedRegion : "all nodes") << " to " << params[3];
					pushEvaluator(params[3], Property::DISPLACEMENT_Y);
					for (size_t e = 0; e < region.size(); e++) {
						region[e]->addSettings(Property::DISPLACEMENT_Y, evaluators.back());
					}
					return;
				}
				if (!params[2].compare(0, 2, "uz")) {
					ESINFO(DETAILS) << "WB: SET uz to " << (_selectedRegion.size() ? "region " + _selectedRegion : "all nodes") << " to " << params[3];
					pushEvaluator(params[3], Property::DISPLACEMENT_Z);
					for (size_t e = 0; e < region.size(); e++) {
						region[e]->addSettings(Property::DISPLACEMENT_Z, evaluators.back());
					}
					return;
				}
				ESINFO(GLOBAL_ERROR) << "Not implemented d option '" << params[2] << "'";
			}
		}
	}

	pushTableEvaluator(evaluators, _tables, params[3]);
	std::vector<Element*> &region = _selectedRegion.size() ? _regions[_selectedRegion] : nodes;
	if (!params[2].compare(0, 2, "ux")) {
		evaluators.back()->property() = Property::DISPLACEMENT_X;
		ESINFO(DETAILS) << "WB: SET ux to " << (_selectedRegion.size() ? "region " + _selectedRegion : "all nodes") << " to table " << evaluators.back()->name();
		for (size_t e = 0; e < region.size(); e++) {
			region[e]->addSettings(Property::DISPLACEMENT_X, evaluators.back());
		}
	}
	if (!params[2].compare(0, 2, "uy")) {
		evaluators.back()->property() = Property::DISPLACEMENT_Y;
		ESINFO(DETAILS) << "WB: SET uy to " << (_selectedRegion.size() ? "region " + _selectedRegion : "all nodes") << " to table " << evaluators.back()->name();
		for (size_t e = 0; e < region.size(); e++) {
			region[e]->addSettings(Property::DISPLACEMENT_Y, evaluators.back());
		}
	}
	if (!params[2].compare(0, 2, "uz")) {
		evaluators.back()->property() = Property::DISPLACEMENT_Z;
		ESINFO(DETAILS) << "WB: SET uz to " << (_selectedRegion.size() ? "region " + _selectedRegion : "all nodes") << " to table " << evaluators.back()->name();
		for (size_t e = 0; e < region.size(); e++) {
			region[e]->addSettings(Property::DISPLACEMENT_Z, evaluators.back());
		}
	}
}

void WorkbenchParser::force(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);

	if (!params[1].compare(0, 3, "all")) {
		if (!params[2].compare(0, 3, "all")) {
			ESINFO(GLOBAL_ERROR) << "Broken WORKBENCH INPUT";
		} else {
			ESINFO(GLOBAL_ERROR) << "Not implemented f option '" << params[2] << "'";
		}
	}

	pushTableEvaluator(evaluators, _tables, params[3]);
	const std::vector<Element*> &region = _regions[params[1]];
	if (!params[2].compare(0, 2, "fx")) {
		evaluators.back()->property() = Property::FORCE_X;
		ESINFO(DETAILS) << "WB: SET fx to region " << params[1] << " to table " << evaluators.back()->name();
		for (size_t e = 0; e < region.size(); e++) {
			region[e]->addSettings(Property::FORCE_X, evaluators.back());
		}
	}
	if (!params[2].compare(0, 2, "fy")) {
		evaluators.back()->property() = Property::FORCE_Y;
		ESINFO(DETAILS) << "WB: SET fy to region " << params[1] << " to table " << evaluators.back()->name();
		for (size_t e = 0; e < region.size(); e++) {
			region[e]->addSettings(Property::FORCE_Y, evaluators.back());
		}
	}
	if (!params[2].compare(0, 2, "fz")) {
		evaluators.back()->property() = Property::FORCE_Z;
		ESINFO(DETAILS) << "WB: SET fz to region " << params[1] << " to table " << evaluators.back()->name();
		for (size_t e = 0; e < region.size(); e++) {
			region[e]->addSettings(Property::FORCE_Z, evaluators.back());
		}
	}
}

void WorkbenchParser::obstacle(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes)
{
	std::vector<std::string> params = divide(_line);

	if (!_selectedRegion.size()) {
		ESINFO(GLOBAL_ERROR) << "You have to set region for obstacle";
	}
	std::vector<Element*> &region = _regions[_selectedRegion];
	evaluators.push_back(new CoordinatesEvaluator(params[1], _mesh.coordinates(), Property::OBSTACLE));
	evaluators.push_back(new CoordinatesEvaluator(params[2], _mesh.coordinates(), Property::NORMAL_DIRECTION));

	ESINFO(DETAILS) << "WB: SET obstacle to region " << _selectedRegion << " to " << params[1] << ", direction: " << params[2];

	for (size_t e = 0; e < region.size(); e++) {
		region[e]->addSettings(Property::OBSTACLE, evaluators[evaluators.size() - 2]);
		region[e]->addSettings(Property::NORMAL_DIRECTION, evaluators[evaluators.size() - 1]);
	}
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
	_selectedRegion = Parser::strip(params[2]);
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
	} else {
		ESINFO(GLOBAL_ERROR) << "Not implemented esel option '" << params[1] << "'";
	}
}




