
#include "options.h"
#include "esbasis.h"

static struct option long_options[] = {
		{"input",  required_argument, 0, 'i'},
		{"path",  required_argument, 0, 'p'},
		{"config",  required_argument, 0, 'c'},
		{"testing", no_argument, 0, 't'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

using namespace espreso;
using namespace input;

Options::Options(int* argc, char*** argv)
: executable((*argv)[0]), verboseLevel(0), testingLevel(0), measureLevel(0)
{
	auto printOption = [] (const std::string &opt, const std::string &desc) {
		if (opt.length() < 16) {
			ESINFO(ALWAYS) << "\t" << opt << "\t\t" << desc;
		} else {
			ESINFO(ALWAYS) << "\t" << opt << "\t" << desc;
		}
	};

	int option_index, option;
	while (true) {
		option = getopt_long(*argc, *argv, "i:p:c:vtmh", long_options, &option_index);
		if (option == -1) {
			break;
		}

		switch (option) {
		case 'v':
			verboseLevel++;
			break;
		case 't':
			testingLevel++;
			break;
		case 'm':
			measureLevel++;
			break;
		case 'i':
			input = std::string(optarg);
			input.erase(0, input.find_first_not_of('='));
			break;
		case 'p':
			path = std::string(optarg);
			path.erase(0, path.find_first_not_of('='));
			break;
		case 'c': {
			std::string config(optarg);
			std::ifstream file(config);
			if (!file.is_open()) {
				ESINFO(ERROR) << "Cannot load configuration file '" << config << "'";
			}
			setFromFile(config);
			break;
		}
		case 'h':
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
			exit(EXIT_SUCCESS);
			break;
		case '?':
			break;
		}
	}

	setFromFile("espreso.config");

	if (!path.size()) {
		if (*argc < 2) { // compatibility with old version of ESPRESO binary
			ESINFO(ERROR) << "Specify path to an example. Run 'espreso -h' for more info.";
		} else {
			path = (*argv)[1];
			verboseLevel = 3;
			measureLevel = 3;
			optind++;
		}
	}

	while (optind < *argc) {
		nameless.push_back(std::string((*argv)[optind++]));
	}

	configure();
}

void Options::setFromFile(const std::string &file)
{
	std::vector<Description> description = {
			{STRING_PARAMETER, "PATH", "A path to an example."},
			{STRING_PARAMETER, "INPUT", "An input format: [generator, matsol, workbench, esdata, openfoam]"}
	};

	Configuration configuration(description, file);

	path = configuration.value("PATH", path);
	input = configuration.value("INPUT", input);
}

static bool caseInsensitiveCmp(char c1, char c2) { return std::tolower(c1) == std::tolower(c2); }

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
	MPI_Comm_rank(MPI_COMM_WORLD, &config::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::MPIsize);
	config::executable = executable;

	std::signal(SIGFPE, signalHandler);
	std::signal(SIGSEGV, signalHandler);

	config::info::verboseLevel += verboseLevel;
	config::info::testingLevel += testingLevel;
	config::info::measureLevel += measureLevel;

	std::vector<std::pair<std::string, config::mesh::Input> > inputs = {
			{ "GENERATOR", config::mesh::GENERATOR },
			{ "MATSOL", config::mesh::ANSYS_MATSOL },
			{ "WORKBENCH", config::mesh::ANSYS_WORKBENCH },
			{ "OPENFOAM", config::mesh::OPENFOAM },
			{ "ESDATA", config::mesh::ESDATA },
	};

	for (size_t i = 0; i < inputs.size(); i++) {
		if (input.size() == inputs[i].first.size() && std::equal(input.begin(), input.end(), inputs[i].first.begin(), caseInsensitiveCmp)) {
			config::mesh::input = inputs[i].second;
			break;
		}
		if (input.size() && i + 1 == inputs.size()) {
			ESINFO(ERROR) << "Unknown input: '" << input << "'";
		}
	}

}

std::ostream& operator<<(std::ostream& os, const Options &options)
{
	os << "input: '" << options.input << "'\n";
	os << "path: '" << options.path << "'\n";
	os << "verbosity level: " << options.verboseLevel << "\n";
	os << "testing level: " << options.testingLevel << "\n";
	os << "time measuring level: " << options.measureLevel << "\n";
	os << "nameless: ";
	for (size_t i = 0; i < options.nameless.size(); i++) {
		os << options.nameless[i] << " ";
	}
	os << "\n";

	return os;
}




