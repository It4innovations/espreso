
#include "options.h"
#include "esbasis.h"

static struct option long_options[] = {
		{"input",  required_argument, 0, 'i'},
		{"path",  required_argument, 0, 'p'},
		{"config",  required_argument, 0, 'c'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

using namespace espreso;
using namespace input;

static void printOption(const std::string &opt, const std::string &desc) {
	if (opt.length() < 16) {
		ESINFO(ALWAYS) << "\t" << opt << "\t\t" << desc;
	} else {
		ESINFO(ALWAYS) << "\t" << opt << "\t" << desc;
	}
};

static void printConfigOption(const Description &option)
{
	auto tabs = [] (size_t size) {
		std::string str;
		for (int i = 0; i < size; i++) {
			str += "\t";
		}
		return str;
	};

	auto defaultValue = [] (const Description &option) {
		std::stringstream ss;
		switch (option.type) {
		case INTEGER_PARAMETER:
			ss << *(static_cast<int*>(option.value));
			break;
		case LONG_PARAMETER:
			ss << *(static_cast<long*>(option.value));
			break;
		case SIZE_PARAMETER:
			ss << *(static_cast<size_t*>(option.value));
			break;
		case DOUBLE_PARAMETER:
			ss << *(static_cast<double*>(option.value));
			break;
		case STRING_PARAMETER:
			ss << *(static_cast<std::string*>(option.value));
			break;
		case BOOLEAN_PARAMETER:
			ss << *(static_cast<bool*>(option.value));
			break;
		}
		return " [" + ss.str() + "] ";
	};

	std::string def = defaultValue(option);
	ESINFO(ALWAYS) << tabs(1) << "  " << option.name << def << tabs((37 - option.name.size() - def.size()) / 8) << "  " << option.description;
	for (size_t i = 0; i < option.options.size(); i++) {
		ESINFO(ALWAYS) << tabs(6) << i << " -> " << option.options[i];
	}
}

static void printFileOptions(size_t level)
{
	switch (level) {
	case 0:
		ESINFO(ALWAYS) << "CONFIGURATION FILE:\n";

		ESINFO(ALWAYS) << "\tThe computation of the ESPRESO depends on various parameters. ";
		ESINFO(ALWAYS) << "\tParameters can be set by a configuration file specified by parameter -c. ";
		ESINFO(ALWAYS) << "\tAll parameters in the file have the following form: ";

		ESINFO(ALWAYS) << "\n\t\tPARAMETER = VALUE\n";

		ESINFO(ALWAYS) << "\tThe ESPRESO accepts the following parameters: [default value]\n";

		for (size_t i = 0; i < config::description.size(); i++) {
			if (config::description[i].writeToHelp == WRITE_TO_HELP) {
				printConfigOption(config::description[i]);
			}
		}
		ESINFO(ALWAYS) << "\n";
		break;
	case 1:
		ESINFO(ALWAYS) << "\t--- Internal parameters ---\n";

		for (size_t i = 0; i < config::description.size(); i++) {
			if (config::description[i].writeToHelp == INGNORE_IN_HELP) {
				printConfigOption(config::description[i]);
			}
		}
		ESINFO(ALWAYS) << "\n";
		break;
	}
}

static bool caseInsensitiveCmp(char c1, char c2) { return std::tolower(c1) == std::tolower(c2); }

Options::Options(int* argc, char*** argv)
: executable((*argv)[0]), verbose(0), measure(0), testing(0), help(0)
{
	int option_index, option;

	std::vector<struct option> opts;
	for (size_t i = 0; i < config::description.size(); i++) {
		opts.push_back({config::description[i].name.c_str(), required_argument, 0, 'd'});
	}

	option_index = 0;
	while (long_options[option_index].name != '\0') {
		opts.push_back(long_options[option_index++]);
	}


	// load configuration first
	while (true) {
		option = getopt_long(*argc, *argv, "i:p:c:vtmh", opts.data(), &option_index);
		if (option == -1) {
			break;
		}

		switch (option) {
		case 'c': {
			std::string config(optarg);
			config::env::configurationFile = config;
			std::ifstream file(config);
			if (!file.is_open()) {
				ESINFO(ERROR) << "Cannot load configuration file '" << config << "'";
			}
			break;
		}
		case 'h':
			if (!help) {
				ESINFO(ALWAYS) << "Usage: espreso [OPTIONS] [PARAMETERS]\n";

				ESINFO(ALWAYS) << "OPTIONS:\n";
				printOption("-h, --help", "show this message");
				printOption("-i, --input=INPUT", "input format: [generator, matsol, workbench, esdata, openfoam]");
				printOption("-p, --path=PATH", "path to an example");
				printOption("-c, --config=FILE", "file with ESPRESO configuration");
				printOption("-v,vv,vvv", "verbose level");
				printOption("-t,tt,ttt", "testing level");
				printOption("-m,mm,mmm", "time measuring level");

				ESINFO(ALWAYS) << "\nPARAMETERS:\n";
				ESINFO(ALWAYS) << "\tlist of nameless parameters for a particular example\n";
			}

			printFileOptions(help++);
			break;
		case '?':
			break;
		}
	}

	if (help) {
		exit(EXIT_SUCCESS);
	}


	Configuration::fill(config::description, config::env::configurationFile);

	optind = 0;
	while (true) {
		option = getopt_long(*argc, *argv, "i:p:c:vtmh", opts.data(), &option_index);
		if (option == -1) {
			break;
		}

		switch (option) {
		case 'd': {
			std::string param(opts[option_index].name);
			for (size_t i = 0; i < config::description.size(); i++) {
				if (param.size() == config::description[i].name.size() && std::equal(param.begin(), param.end(), config::description[i].name.begin(), caseInsensitiveCmp)) {
					switch (config::description[i].type) {
					case INTEGER_PARAMETER:
						*(static_cast<int*>(config::description[i].value)) = std::stoi(optarg);
						break;
					case LONG_PARAMETER:
						*(static_cast<long*>(config::description[i].value)) = std::stol(optarg);
						break;
					case SIZE_PARAMETER:
						*(static_cast<size_t*>(config::description[i].value)) = std::stoull(optarg);
						break;
					case DOUBLE_PARAMETER:
						*(static_cast<double*>(config::description[i].value)) = std::stod(optarg);
						break;
					case STRING_PARAMETER:
						*(static_cast<std::string*>(config::description[i].value)) = optarg;
						break;
					case BOOLEAN_PARAMETER:
						*(static_cast<bool*>(config::description[i].value)) = std::stoi(optarg);
						break;
					}
					break;
				}
			}
			break;
		}
		case 'v':
			verbose++;
			break;
		case 't':
			testing++;
			break;
		case 'm':
			measure++;
			break;
		case 'i':
			input = std::string(optarg);
			input.erase(0, input.find_first_not_of('='));
			break;
		case 'p':
			path = std::string(optarg);
			path.erase(0, path.find_first_not_of('='));
			break;
		case 'c':
			// configuration is already loaded
			break;
		case '?':
			break;
		}
	}

	config::info::verboseLevel += verbose;
	config::info::measureLevel += measure;
	config::info::testingLevel += testing;

	if (path.size()) {
		config::mesh::PATH = path;
	} else {
		path = config::mesh::PATH;
	}

	if (!path.size()) {
		if (*argc < 2) { // compatibility with old version of ESPRESO binary
			ESINFO(ERROR) << "Specify path to an example. Run 'espreso -h' for more info.";
		} else {
			path = (*argv)[1];
			config::info::verboseLevel = 3;
			config::info::measureLevel = 3;
			optind++;
		}
	}

	while (optind < *argc) {
		nameless.push_back(std::string((*argv)[optind++]));
	}

	configure();
}

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGSEGV:
		ESINFO(ERROR) << "Invalid memory reference";
		break;
	case SIGFPE:
		ESINFO(ERROR) << "Erroneous arithmetic operation";
		break;
	}
}

void Options::configure()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::env::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::env::MPIsize);
	config::env::executable = executable;

	std::signal(SIGFPE, signalHandler);
	std::signal(SIGSEGV, signalHandler);

	std::vector<std::pair<std::string, config::mesh::INPUTalternatives> > inputs = {
			{ "GENERATOR", config::mesh::GENERATOR },
			{ "MATSOL", config::mesh::MATSOL },
			{ "WORKBENCH", config::mesh::WORKBENCH },
			{ "OPENFOAM", config::mesh::OPENFOAM },
			{ "ESDATA", config::mesh::ESDATA },
	};

	for (size_t i = 0; i < inputs.size(); i++) {
		if (input.size() == inputs[i].first.size() && std::equal(input.begin(), input.end(), inputs[i].first.begin(), caseInsensitiveCmp)) {
			config::mesh::INPUT = inputs[i].second;
			break;
		}
		if (input.size() && i + 1 == inputs.size()) {
			ESINFO(ERROR) << "Unknown input: '" << input << "'";
		}
	}

}

std::ostream& espreso::operator<<(std::ostream& os, const Options &options)
{
	os << "input: '" << options.input << "'\n";
	os << "path: '" << options.path << "'\n";
	os << "nameless: ";
	for (size_t i = 0; i < options.nameless.size(); i++) {
		os << options.nameless[i] << " ";
	}
	os << "\n";

	return os;
}




