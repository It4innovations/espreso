
#include "logging.h"

namespace eslog {

std::vector<Checkpoint> Log::checkpoints = { Checkpoint("Start", Log::time(), 0) };

Log::Log(LogEvent event): event(event)
{
	auto indent = [&] (int tabs) { for (int t = 0; t < tabs; t++) { os << "  "; } };

	switch (event) {
	case CHECKPOINT3:
		indent(2);
		break;
	case CHECKPOINT2:
		indent(1);
		break;
	case ERROR:
		os << "ESPRESO ERROR : ";
		break;
	}
}

Log::~Log()
{
	if (event == ERROR) {
		fprintf(stderr, "%s", os.str().c_str());
		fprintf(stderr, "ESPRESO EXITED WITH ERROR ON PROCESS %d.\n", esconfig::MPIrank);
		fflush(stderr);
		exit(EXIT_FAILURE);
	}

	switch (event) {
	case CHECKPOINT3:
		checkpoints.push_back(Checkpoint(os.str(), time(), 3));
		break;
	case CHECKPOINT2:
		checkpoints.push_back(Checkpoint(os.str(), time(), 2));
		break;
	case CHECKPOINT1:
		checkpoints.push_back(Checkpoint(os.str(), time(), 1));
		break;
	case SUMMARY:
		checkpoints.push_back(Checkpoint(os.str(), time(), 0));
		evaluateCheckpoints();
		return;
	}

	os << std::endl;

	if (esconfig::MPIrank != 0) {
		return; // only first process print results
	}

	fprintf(stdout, "%s", os.str().c_str());
	fflush(stdout);
}

void Log::evaluateCheckpoints()
{
	if (esconfig::MPIrank != 0) {
		return;
	}

	for (size_t i = 0; i < checkpoints.size(); i++) {
		std::cout << std::fixed << std::setprecision(3) << checkpoints[i].name << " : " << checkpoints[i].time << "\n";
	}
}

Test::~Test()
{
	os << " : " << (error ? "FAILED" : "PASSED") << std::endl;

	if (error) {
		ESLOG(ERROR) << os.str();
	}

	fprintf(stdout, "%s", os.str().c_str());
	fflush(stdout);
}
}



