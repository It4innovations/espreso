
#include "logging.h"

namespace espreso {

std::vector<Checkpoint> Measure::checkpoints = { Checkpoint("Start", Measure::time(), 0) };

Test::~Test()
{
	os << " : " << (error ? "FAILED" : "PASSED") << std::endl;

	if (error) {
		ESINFO(ERROR) << os.str();
	}

	fprintf(stdout, "%s", os.str().c_str());
	fflush(stdout);
}


Info::~Info()
{
	if (_plain) {
		if (config::MPIrank == 0) {
			std::cout << os.str();
		}
		return;
	}
	if (event == ERROR) {
		fprintf(stderr, "%s\n", os.str().c_str());
		fprintf(stderr, "ESPRESO EXITED WITH ERROR ON PROCESS %d.\n", config::MPIrank);
		fflush(stderr);
		exit(EXIT_FAILURE);
	}

	os << std::endl;

	if (config::MPIrank != 0) {
		return; // only first process print results
	}

	fprintf(stdout, "%s", os.str().c_str());
	fflush(stdout);
}


Measure::~Measure()
{
//	switch (event) {
//	case CHECKPOINT3:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 3));
//		break;
//	case CHECKPOINT2:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 2));
//		break;
//	case CHECKPOINT1:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 1));
//		break;
//	case SUMMARY:
//		checkpoints.push_back(Checkpoint(os.str(), time(), 0));
//		evaluateCheckpoints();
//		return;
//	}

	os << std::endl;

	if (config::MPIrank != 0) {
		return; // only first process print results
	}

	fprintf(stdout, "%s", os.str().c_str());
	fflush(stdout);
}

void Measure::evaluateCheckpoints()
{
	if (config::MPIrank != 0) {
		return;
	}

	for (size_t i = 0; i < checkpoints.size(); i++) {
		std::cout << std::fixed << std::setprecision(3) << checkpoints[i].name << " : " << checkpoints[i].time << "\n";
	}
}

double Measure::processMemory()
{
	std::ifstream file("/proc/self/status");
	eslocal result = -1;
	std::string line, label("VmRSS:");

	while (getline(file, line)) {
		if (line.find(label.c_str(), 0, label.size()) == 0) {
			file.close();
			std::stringstream(line) >> label >> result;
			return result / 1024.0;
		}
	}
	file.close();
	return 0;
}

double Measure::usedRAM()
{
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	return ((memInfo.totalram - memInfo.freeram) * memInfo.mem_unit) / 1024.0 / 1024.0;
}

double Measure::availableRAM()
{
	struct sysinfo memInfo;
	sysinfo(&memInfo);
	return (memInfo.totalram * memInfo.mem_unit) / 1024.0 / 1024.0;
}

}



