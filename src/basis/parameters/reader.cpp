
#include "reader.h"

using namespace espreso;

static struct option long_options[] = {
		{"input",  required_argument, 0, 'i'},
		{"path",  required_argument, 0, 'p'},
		{"config",  required_argument, 0, 'c'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

ParametersReader::ParametersReader(const std::vector<Parameter> &parameters)
: _parameters(parameters)
{
	std::sort(_parameters.begin(), _parameters.end());
}

ArgsConfiguration ParametersReader::fromArguments(int *argc, char*** argv, const std::vector<Parameter> &params)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);

	ParametersReader reader(params);
	ArgsConfiguration conf;

	int option_index, option;
	std::string options("i:p:c:vtmh");

	std::vector<struct option> opts;
	for (size_t i = 0; i < params.size(); i++) {
		opts.push_back({params[i].name.c_str(), required_argument, 0, 'd'});
	}
	option_index = 0;
	while (long_options[option_index].name != '\0') {
		opts.push_back(long_options[option_index++]);
	}

	// check for help
	size_t helpVerboseLevel = 0;
	while ((option = getopt_long(*argc, *argv, options.c_str(), opts.data(), &option_index)) != -1) {
		if (option == 'h') {
			helpVerboseLevel++;
		}
	}
	if (helpVerboseLevel) {
		printHelp(helpVerboseLevel + 1);
		exit(0);
	}

	// read nameless parameters
	while (optind < *argc) {
		conf.nameless.push_back(std::string((*argv)[optind++]));
	}

	// load configuration before other parameters
	optind = 0;
	bool configured = false;
	while ((option = getopt_long(*argc, *argv, options.c_str(), opts.data(), &option_index)) != -1) {
		if (option == 'c') {
			conf.path = optarg;
			conf = reader.read(conf, 1);
			configured = true;
		}
	}
	if (!configured && std::ifstream(config::env::configurationFile).good()) {
		// read the default configuration file
		conf.path = config::env::configurationFile;
		conf = reader.read(conf, 1);
	}

	// read the rest parameters
	optind = 0;
	while ((option = getopt_long(*argc, *argv, "i:p:c:vtmh",opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'd':
			reader._setParameter(opts[option_index].name, optarg);
			break;
		case 'v':
			config::info::VERBOSE_LEVEL++;
			break;
		case 't':
			config::info::TESTING_LEVEL++;
			break;
		case 'm':
			config::info::MEASURE_LEVEL++;
			break;
		case 'i':
			reader._setParameter("INPUT", optarg);
			break;
		case 'p':
			reader._setParameter("PATH", optarg);
			break;
		}
	}

	// compatibility with old version of ESPRESO
	conf.path = config::mesh::PATH;
	if (conf.path.size() == 0 && conf.nameless.size() + 1 == (size_t)(*argc)) {
		if (conf.nameless.size()) {
			conf.path = conf.nameless.front();
			conf.nameless.erase(conf.nameless.begin());
			config::info::VERBOSE_LEVEL = 3;
			config::info::MEASURE_LEVEL = 3;
		}
	}
	if (!conf.path.size()) {
		ESINFO(GLOBAL_ERROR) << "Specify path to an example. Run 'espreso -h' for more info.";
	}

	return conf;
}

ArgsConfiguration ParametersReader::fromConfigurationFile(const ArgsConfiguration &conf, const std::vector<Parameter> &params)
{
	ParametersReader reader(params);
	return reader.read(conf, 1);
}

ArgsConfiguration ParametersReader::fromConfigurationFileWOcheck(const ArgsConfiguration &conf, const std::vector<Parameter> &params)
{
	ParametersReader reader(params);
	return reader.read(conf, 0);
}


ArgsConfiguration ParametersReader::read(const ArgsConfiguration &configuration, size_t verboseLevel)
{
	std::ifstream file(configuration.path);
	ArgsConfiguration conf(configuration);

	std::vector<std::string> cmdArgs, attributes, values;

	if (!file.is_open()) {
		ESINFO(ERROR) << "A configuration file on path '" << configuration.path << "' not found.";
	}

	std::string line;
	while (getline(file, line, '\n')) {
		std::string parameter = Parser::getParameter(line);
		if (parameter.size() == 0) {
			// skip empty line
			continue;
		}
		if (Parser::getValue(line) == "{") {
			attributes.clear();
			values.clear();
			while (getline(file, line, '\n')) {
				std::string attribute = Parser::getParameter(line);
				if (attribute.size() == 0) {
					// skip empty line
					continue;
				}
				if (attribute == "}") {
					if (!_setParameters(parameter, attributes, values) && verboseLevel) {
						ESINFO(ALWAYS) << Info::TextColor::RED << "Unknown parameter '" << parameter << "'";
					}
					break;
				} else {
					attributes.push_back(Parser::getParameter(line));
					values.push_back(Parser::getValue(line));
				}
			}
			continue;
		}

		if (StringCompare::caseInsensitiveEq(parameter, "CMD_LINE_ARGUMENTS")) {
			cmdArgs = Parser::split(Parser::getValue(line), " ");
			if (conf.nameless.size() < cmdArgs.size()) {
				ESINFO(GLOBAL_ERROR) << "Too few nameless command line arguments. \n";
			}
		} else {
			if (!_setParameter(parameter, Parser::getValue(line)) && verboseLevel) {
				ESINFO(ALWAYS) << Info::TextColor::RED << "Unknown parameter '" << parameter << "'";
			}
		}
	}

	for (size_t i = 0; i < cmdArgs.size(); i++) {
		if (!_setParameter(cmdArgs[i], conf.nameless[i]) && verboseLevel) {
			ESINFO(ALWAYS) << "Unknown command line parameter '" << cmdArgs[i] << "'";
		}
	}
	// remove already red nameless parameters
	conf.nameless.erase(conf.nameless.begin(), conf.nameless.begin() + cmdArgs.size());

	file.close();
	return conf;
}

bool ParametersReader::_setParameter(const std::string &parameter, const std::string &value)
{
	auto it = std::lower_bound(_parameters.begin(), _parameters.end(), parameter);
	if (it != _parameters.end() && StringCompare::caseInsensitiveEq(it->name, parameter)) {
		it->set(value);
		return true;
	} else {
		return false;
	}
}

bool ParametersReader::_setParameters(const std::string &parameter, const std::vector<std::string> &attributes, const std::vector<std::string> &values)
{
	auto it = std::lower_bound(_parameters.begin(), _parameters.end(), parameter);
	if (it != _parameters.end() && StringCompare::caseInsensitiveEq(it->name, parameter)) {
		for (size_t i = 0; i < attributes.size(); i++) {
			it->set(attributes[i], values[i]);
		}
		return true;
	} else {
		return false;
	}
}

static void printParameter(const Parameter &parameter)
{
	auto tabs = [] (int size) {
		std::string str;
		for (int i = 0; i < size; i++) {
			str += "\t";
		}
		return str;
	};

	ESINFO(ALWAYS)
			<< tabs(1) << "  " << parameter.name
			<< " [" << parameter.get() << "]"
			<< tabs((34 - parameter.name.size() - parameter.get().size()) / 8)
			<< "  " << parameter.description;

	for (size_t i = 0; i < parameter.data->options(); i++) {
		ESINFO(ALWAYS) << tabs(6) << "[" << i << "] " << parameter.data->optionName(i) << " - " << parameter.data->optionDesc(i);
	}
}

void ParametersReader::printHelp(size_t verboseLevel)
{
	auto printOption = [](const std::string &opt, const std::string &desc) {
		if (opt.length() < 16) {
			ESINFO(ALWAYS) << "\t" << opt << "\t\t" << desc;
		} else {
			ESINFO(ALWAYS) << "\t" << opt << "\t" << desc;
		}
	};

	ESINFO(ALWAYS) << "Usage: espreso [OPTIONS] [PARAMETERS]\n";

	ESINFO(ALWAYS) << "OPTIONS:\n";
	printOption("-h, --help", "show this message");
	printOption("-i, --input=INPUT", "input format: [MATSOL, WORKBENCH, OPENFOAM, ESDATA, GENERATOR]");
	printOption("-p, --path=PATH", "path to an example");
	printOption("-c, --config=FILE", "file with ESPRESO configuration");
	printOption("-v,vv,vvv", "verbose level");
	printOption("-t,tt,ttt", "testing level");
	printOption("-m,mm,mmm", "time measuring level");

	printOption("--PARAMETER=VALUE", "rewrite an arbitrary parameter from a configuration file");

	ESINFO(ALWAYS) << "\nPARAMETERS:\n";
	ESINFO(ALWAYS) << "\tlist of nameless parameters for a particular example\n";

	ESINFO(ALWAYS) << "CONFIGURATION FILE PARAMETERS:\n";

	ESINFO(ALWAYS) << "\tThe computation of the ESPRESO depends on various parameters. ";
	ESINFO(ALWAYS) << "\tParameters can be set by a configuration file specified by parameter -c. ";
	ESINFO(ALWAYS) << "\tAll parameters in the file have the following form: ";

	ESINFO(ALWAYS) << "\n\t\tPARAMETER = VALUE\n";

	ESINFO(ALWAYS) << "\tThe ESPRESO accepts the following parameters: [default value]\n";

	printParametersHelp(config::parameters, verboseLevel);
}

void ParametersReader::printParametersHelp(const std::vector<Parameter> &params, size_t verboseLevel)
{
	for (size_t i = 0; i < verboseLevel; i++) {
		size_t n = 0;
		for (size_t j = 0; j < params.size(); j++) {
			if (params[j].verboseLevel == i) {
				printParameter(params[j]);
				n++;
			}
		}
		if (n) {
			ESINFO(ALWAYS);
		}
	}
}

void ParametersReader::printParameters(const std::vector<Parameter> &params, size_t verboseLevel)
{
	for (size_t i = 0; i < verboseLevel; i++) {
		size_t n = 0;
		for (size_t j = 0; j < params.size(); j++) {
			if (params[j].verboseLevel == i) {
				auto attributes = params[j].data->attributes();
				if (attributes.size()) {
					ESINFO(ALWAYS) << "\t" << params[j].name << " == {";
					for (size_t k = 0; k < attributes.size(); k++) {
						ESINFO(ALWAYS) << "\t  " << attributes[k] << " == " << params[j].get(attributes[k]);
					}
					ESINFO(ALWAYS) << "\t}";
				} else {
					ESINFO(ALWAYS) << "\t" << params[j].name << " == " << params[j].get();
				}
				n++;
			}
		}
		if (n) {
			ESINFO(ALWAYS);
		}
	}
}


