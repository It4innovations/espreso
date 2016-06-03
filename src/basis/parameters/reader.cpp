
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

Configuration ParametersReader::arguments(int *argc, char*** argv, const std::vector<Parameter> &params)
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);

	ParametersReader reader(params);
	Configuration conf;

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
	int help_level = 0;
	while ((option = getopt_long(*argc, *argv, options.c_str(), opts.data(), &option_index)) != -1) {
		if (option == 'h') {
			help_level++;
		}
	}
	if (help_level) {
		printHelp(help_level);
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
			conf = reader.read(conf);
			configured = true;
		}
	}
	if (!configured) {
		// read the default configuration file
		conf.path = config::env::configurationFile;
		conf = reader.read(conf);
	}

	// read the rest parameters
	optind = 0;
	while ((option = getopt_long(*argc, *argv, "i:p:c:vtmh",opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'd':
			reader.setParameter(opts[option_index].name, optarg);
			break;
		case 'v':
			config::info::verboseLevel++;
			break;
		case 't':
			config::info::testingLevel++;
			break;
		case 'm':
			config::info::measureLevel++;
			break;
		case 'i':
			reader.setParameter("INPUT", optarg);
			break;
		case 'p':
			reader.setParameter("PATH", optarg);
			break;
		}
	}

	// compatibility with old version of ESPRESO
	conf.path = config::mesh::PATH;
	if (conf.path.size() == 0 && conf.nameless.size() + 1 == *argc) {
		if (conf.nameless.size()) {
			conf.path = conf.nameless.front();
			conf.nameless.erase(conf.nameless.begin());
		} else {
			ESINFO(GLOBAL_ERROR) << "Specify path to an example. Run 'espreso -h' for more info.";
		}
	}

	return conf;
}

Configuration ParametersReader::configuration(const Configuration &conf, const std::vector<Parameter> &params)
{
	ParametersReader reader(params);
	return reader.read(conf);
}

void ParametersReader::pickParameter(const Configuration &conf, const std::string parameter, const std::vector<Parameter> &params)
{
	ParametersReader reader(params);

	std::ifstream file(conf.path);
	if (!file.is_open()) {
		ESINFO(ERROR) << "A configuration file on path '" << conf.path << "' not found.";
	}

	std::string line;
	while (getline(file, line, '\n')) {
		std::string p = Parser::getParameter(line);
		if (StringCompare::caseInsensitiveEq(p, parameter)) {
			reader.setParameter(p, Parser::getValue(line));
			return;
		}
	}
}


Configuration ParametersReader::read(const Configuration &configuration)
{
	std::ifstream file(configuration.path);
	Configuration conf(configuration);

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

		if (StringCompare::caseInsensitiveEq(parameter, "CMD_LINE_ARGUMENTS")) {
			std::vector<std::string> params = Parser::split(Parser::getValue(line), " ");
			if (conf.nameless.size() < params.size()) {
				ESINFO(GLOBAL_ERROR) << "Too few nameless command line arguments. \n";
			}
			for (size_t i = 0; i < params.size(); i++) {
				setParameter(params[i], conf.nameless[i]);
			}
			// remove already red nameless parameters
			conf.nameless.erase(conf.nameless.begin(), conf.nameless.begin() + params.size());
		} else {
			setParameter(parameter, Parser::getValue(line));
		}
	}

	file.close();
	return conf;
}

void ParametersReader::setParameter(const std::string &parameter, const std::string &value)
{
	auto it = std::lower_bound(_parameters.begin(), _parameters.end(), parameter);
	if (it != _parameters.end() && StringCompare::caseInsensitiveEq(it->name, parameter)) {
		it->set(value);
	} else {
		ESINFO(ALWAYS) << "Unknown parameter '" << parameter << "'";
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

void ParametersReader::printHelp(int level)
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

	printParameters(config::parameters, level);
}

void ParametersReader::printParameters(const std::vector<Parameter> &params, int level)
{
	if (level < 1) {
		return;
	}

	for (size_t i = 0; i < params.size(); i++) {
		if (params[i].help == Parameter::Help::WRITE) {
			printParameter(params[i]);
		}
	}
	ESINFO(ALWAYS) << "\n";

	if (level < 2) {
		return;
	}
	ESINFO(ALWAYS) << "\t--- Internal parameters ---\n";

	for (size_t i = 0; i < params.size(); i++) {
		if (params[i].help == Parameter::Help::INGNORE) {
			printParameter(params[i]);
		}
	}
	ESINFO(ALWAYS) << "\n";
}


