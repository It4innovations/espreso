
#include "logging.h"

namespace eslog {

double Log::start = TimeEvent::time();
std::vector<double> Log::lastTimes = { Log::start };

Log::Log(Event event): event(event), error(false)
{
	os << std::fixed << std::setprecision(8);

	auto indent = [&] (int tabs) { for (int t = 0; t < tabs; t++) { os << "  "; } };
	auto time = [&] (int level) {
		if (lastTimes.size() > level) { lastTimes.resize(level); }
		for (int l = lastTimes.size(); l < level; l++) {
			lastTimes.push_back(lastTimes.back());
		}
	};

	switch (event) {
	case CHECKPOINT3:
		time(3);
		indent(2);
		os << (esconfig::info::testingMode ? "CHECKPOINT : " : "");
		break;
	case CHECKPOINT2:
		time(2);
		indent(1);
		os << (esconfig::info::testingMode ? "CHECKPOINT : " : "");
		break;
	case CHECKPOINT1:
		time(1);
		os << (esconfig::info::testingMode ? "CHECKPOINT : " : "");
		break;
	case DURATION:
		os << (esconfig::info::testingMode ? "DURATION : " : "");
		break;
	case TEST_SIMPLE: case TEST_EXPENSIVE:
		os << (esconfig::info::testingMode ? "TEST : " : "");
		break;
	case ERROR:
		os << "ESPRESO ERROR : ";
		error = true;
		break;
	}
}

static void duration(std::ostringstream &os, std::vector<double> lastTimes, int level, bool collective)
{
	double last = lastTimes[level - 1];
	lastTimes[level - 1] = TimeEvent::time();
	double diff = lastTimes[level - 1] - last;
	if (collective) {
		std::vector<double> times(esconfig::MPIsize);
		MPI_Gather(&diff, sizeof(double), MPI_BYTE, times.data(), sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

		auto result = std::minmax_element(times.begin(), times.end());;
		os << *result.first << " : " << *result.second;
	} else {
		os << diff;
	}
}

Log::~Log()
{
	if (event != ERROR) {
		os << " : ";
	}

	switch (event) {
	case CHECKPOINT3:
		duration(os ,lastTimes, 3, esconfig::info::verboseLevel > 1);
		break;
	case CHECKPOINT2:
		duration(os ,lastTimes, 2, esconfig::info::verboseLevel > 1);
		break;
	case CHECKPOINT1:
		duration(os ,lastTimes, 1, esconfig::info::verboseLevel > 1);
		break;
	case TEST_SIMPLE: case TEST_EXPENSIVE:
		os << (error ? "FAILED" : "PASSED");
		break;
	case DURATION:
		if (esconfig::info::verboseLevel > 1) {
			MPI_Barrier(MPI_COMM_WORLD);
		}
		os << TimeEvent::time() - start;
	}

	os << std::endl;

	if (esconfig::MPIrank != 0 && !error) {
		return; // only first process print results
	}

	if (error) {
		fprintf(stderr, "%s", os.str().c_str());
		fprintf(stderr, "ESPRESO EXITED WITH ERROR ON PROCESS %d.\n", esconfig::MPIrank);
		fflush(stderr);
		exit(EXIT_FAILURE);
	} else {
		fprintf(stdout, "%s", os.str().c_str());
		fflush(stdout);
	}
}

}



