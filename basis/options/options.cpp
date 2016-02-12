
#include "options.h"

static struct option long_options[] = {
		{"input",  required_argument, 0, 'i'},
		{"path",  required_argument, 0, 'p'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

Options::Options(int* argc, char*** argv): verboseLevel(VERBOSE)
{
	auto printOption = [] (const std::string &opt, const std::string &desc) {
		std::cout << "\t" << opt;
		if (opt.length() < 16) {
			std::cout << "\t\t" << desc << "\n";
		} else {
			std::cout << "\t" << desc << "\n";
		}
	};

	int option_index, option;
	while (true) {
		option = getopt_long(*argc, *argv, "i:p:vh", long_options, &option_index);
		if (option == -1) {
			break;
		}

		switch (option) {
		case 'v':
			verboseLevel++;
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
			std::cout << "Usage: espreso [OPTIONS] [PARAMETERS]\n";

			std::cout << "\nOPTIONS:\n";
			printOption("-h, --help", "show this message");
			printOption("-i, --input=INPUT", "input format: [generator, matsol, workbench, esdata]");
			printOption("-p, --path=PATH", "path to an example");
			printOption("-v,vv,vvv", "verbose level");

			std::cout << "\nPARAMETERS:\n";
			std::cout << "\tlist of nameless parameters for a particular example\n";
			exit(EXIT_SUCCESS);
			break;
		case '?':
			break;
		}
	}

	if (optind == 1 || (optind > 1 && !path.size())) { // compatibility with old version of ESPRESO binary
		if (*argc < 2) {
			std::cerr << "ESPRESO Error: specify path to an example. Run 'espreso -h' for more info.\n";
			exit(EXIT_FAILURE);
		}
		path = (*argv)[1];
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
	MPI_Comm_rank(MPI_COMM_WORLD, &esconfig::MPIrank);
	MPI_Comm_size(MPI_COMM_WORLD, &esconfig::MPIsize);

	esconfig::info::verboseLevel = verboseLevel;

	std::vector<std::pair<std::string, esconfig::mesh::Input> > inputs = {
			{ "GENERATOR", esconfig::mesh::GENERATOR },
			{ "MATSOL", esconfig::mesh::ANSYS_MATSOL },
			{ "WORKBENCH", esconfig::mesh::ANSYS_WORKBENCH },
			{ "OPENFOAM", esconfig::mesh::OPENFOAM },
			{ "ESDATA", esconfig::mesh::ESDATA_IN },
	};
	for (size_t i = 0; i < inputs.size(); i++) {
		if (!std::lexicographical_compare(
				input.begin(), input.end(),
				inputs[i].first.begin(), inputs[i].first.end(),
				caseInsensitiveCmp)) {

			esconfig::mesh::input = inputs[i].second;
		}
	}

}

std::ostream& operator<<(std::ostream& os, const Options &options)
{
	os << "input: '" << options.input << "'\n";
	os << "path: '" << options.path << "'\n";
	os << "verbosity level: " << options.verboseLevel << "\n";
	os << "nameless: ";
	for (size_t i = 0; i < options.nameless.size(); i++) {
		os << options.nameless[i] << " ";
	}
	os << "\n";

	return os;
}




