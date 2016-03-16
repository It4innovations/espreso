
#include "options.h"

static struct option long_options[] = {
		{"input",  required_argument, 0, 'i'},
		{"path",  required_argument, 0, 'p'},
		{"testing", no_argument, 0, 't'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

using namespace espreso;

Options::Options(int* argc, char*** argv)
: verboseLevel(0), testingLevel(0), measureLevel(0)
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
		option = getopt_long(*argc, *argv, "i:p:vtmh", long_options, &option_index);
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
		case 'h':
			ESINFO(ALWAYS) << "Usage: espreso [OPTIONS] [PARAMETERS]\n";

			ESINFO(ALWAYS) << "\nOPTIONS:\n";
			printOption("-h, --help", "show this message");
			printOption("-i, --input=INPUT", "input format: [generator, matsol, workbench, esdata, openfoam]");
			printOption("-p, --path=PATH", "path to an example");
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

	if (optind == 1 || (optind > 1 && !path.size())) { // compatibility with old version of ESPRESO binary
		if (*argc < 2) {
			ESINFO(ERROR) << "Specify path to an example. Run 'espreso -h' for more info.";
		}
		path = (*argv)[1];
		verboseLevel = 3;
		measureLevel = 3;
		optind++;
	}

	while (optind < *argc) {
		nameless.push_back(std::string((*argv)[optind++]));
	}

	configure();
}

static bool caseInsensitiveCmp(char c1, char c2) { return std::tolower(c1) < std::tolower(c2); }

void Options::configure()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &config::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &config::MPIsize);

	config::info::verboseLevel = verboseLevel;
	config::info::testingLevel = testingLevel;
	config::info::measureLevel = measureLevel;

	std::vector<std::pair<std::string, config::mesh::Input> > inputs = {
			{ "GENERATOR", config::mesh::GENERATOR },
			{ "MATSOL", config::mesh::ANSYS_MATSOL },
			{ "WORKBENCH", config::mesh::ANSYS_WORKBENCH },
			{ "OPENFOAM", config::mesh::OPENFOAM },
			{ "ESDATA", config::mesh::ESDATA },
	};
	for (size_t i = 0; i < inputs.size(); i++) {
		if (!std::lexicographical_compare(
				input.begin(), input.end(),
				inputs[i].first.begin(), inputs[i].first.end(),
				caseInsensitiveCmp)) {

			config::mesh::input = inputs[i].second;
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




