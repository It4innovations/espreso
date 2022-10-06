
#include "systeminfo.h"
#include "eslog.hpp"
#include "ecfinfo.h"
#include "envinfo.h"
#include "wrappers/mpi/communication.h"

#include <csignal>
#include <cstdlib>
#include <cstring>
#include <thread>
#include <fstream>
#include <sched.h>
#include <vector>

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGTERM:
		espreso::eslog::failure("SIGTERM - termination request, sent to the program.\n");
		break;
	case SIGSEGV:
		espreso::eslog::failure("SIGSEGV - invalid memory access (segmentation fault).\n");
		break;
	case SIGINT:
		espreso::eslog::failure("SIGINT - external interrupt, usually initiated by the user.\n");
		break;
	case SIGILL:
		espreso::eslog::failure("SIGILL - invalid program image, such as invalid instruction.\n");
		break;
	case SIGABRT:
		espreso::eslog::failure("SIGABRT - abnormal termination condition, as is e.g. initiated by std::abort().\n");
		break;
	case SIGFPE:
		espreso::eslog::failure("SIGFPE - erroneous arithmetic operation such as divide by zero.\n");
		break;
	default:
		espreso::eslog::failure("ESPRESO trigger an error.\n");
	}
	exit(EXIT_FAILURE);
}

namespace espreso {
namespace info {
namespace system {

OPERATIONSYSTEM os()
{
	return OPERATIONSYSTEM::UNIX;
}

const char* commit()
{
#ifdef __ESCOMMIT__
	return __ESCOMMIT__;
#else
	return "UNKNOWN";
#endif
}

const char* cxx()
{
#ifdef __ESCXX__
	return __ESCXX__;
#else
	return "UNKNOWN";
#endif
}

const char* buildpath()
{
#ifdef __ESBUILDPATH__
	return __ESBUILDPATH__;
#else
	return ".";
#endif
}

const char* cxxflags()
{
#ifdef __ESCXXFLAGS__
	return __ESCXXFLAGS__;
#else
	return "UNKNOWN";
#endif
}

int MPI_Init_thread_provided_level()
{
#ifdef __MPI_Init_thread_level__
	return __MPI_Init_thread_level__;
#else
	return MPI_THREAD_SINGLE;
#endif
}

int hwthreads()
{
	return std::thread::hardware_concurrency();
}

CPUInfo cpuinfo()
{
	CPUInfo info;

	std::ifstream is("/proc/cpuinfo");

	int siblings = 0;
	std::string param, colon, value;
	while (is.good()) {
		getline(is, param, '\t');
		getline(is, colon, ':');
		getline(is, value, '\n');
		if (strcmp(param.c_str(), "model name") == 0) {
			strcpy(info.modelName, value.c_str());
		}
		if (strcmp(param.c_str(), "physical id") == 0) {
			info.sockets = std::max(std::stoi(value.c_str(), nullptr, 10), info.sockets);
		}
		if (strcmp(param.c_str(), "siblings") == 0) {
			siblings = std::stoi(value.c_str(), nullptr, 10);
		}
		if (strcmp(param.c_str(), "cpu cores") == 0) {
			info.cores = std::stoi(value.c_str(), nullptr, 10);
		}
	}
	info.hyperthreading = siblings != info.cores;
	info.sockets++;

	return info;
}

long pinnedDomainSize()
{
	cpu_set_t affinity;
	sched_getaffinity(0, sizeof(cpu_set_t), &affinity);
	return CPU_COUNT(&affinity);
}

bool pinningIntersection()
{
	cpu_set_t affinity;
	sched_getaffinity(0, sizeof(cpu_set_t), &affinity);

	std::vector<int> count(hwthreads());
	for (size_t i = 0; i < count.size(); ++i) {
		count[i] = CPU_ISSET(i, &affinity) ? 1 : 0;
	}
	Communication::allReduce(count.data(), nullptr, count.size(), MPI_INT, MPI_SUM, &MPITools::node->within);
	for (size_t i = 0; i < count.size(); ++i) {
		if (count[i] > 1) {
			return true;
		}
	}
	return false;
}

long MemTotal, MemFree, MemAvailable;

long memoryTotal()
{
	if (MemTotal == 0) { // total memory is persistent
		std::ifstream meminfo("/proc/meminfo");
		std::string name, unit;

		meminfo >> name >> MemTotal >> unit;
		meminfo >> name >> MemFree >> unit;
		meminfo >> name >> MemAvailable >> unit;
	}
	return MemTotal;
}

long memoryAvail()
{
	std::ifstream meminfo("/proc/meminfo");
	std::string name, unit;

	meminfo >> name >> MemTotal >> unit;
	meminfo >> name >> MemFree >> unit;
	meminfo >> name >> MemAvailable >> unit;
	return MemAvailable;
}

const char* simd()
{
#if defined(__AVX512F__) && defined(__AVX512DQ__)
	return "AVX-512";
#elif defined(__AVX__)
	return "AVX";
#elif defined(__SSE2__)
	return "SSE2";
#elif defined(__SSE__)
	return "SSE2";
#else
	return "UNKNOWN";
#endif
}

void print()
{
	int ppn = MPITools::node->within.size;
	int threads = info::system::hwthreads();
	int nodes = info::mpi::size / ppn;
	const char* simd = info::system::simd();
	auto cpuinfo = info::system::cpuinfo();
	int secondRankNode = info::mpi::rank == 1 ? MPITools::node->across.rank : 0;
	Communication::allReduce(&secondRankNode, nullptr, 1, MPI_INT, MPI_MAX);

	double memTotal = 0, memAvail = 0;
	if (info::mpi::rank == 0) {
		memAvail = info::system::memoryAvail();
		memTotal = info::system::memoryTotal();
	}
	int asynchronousThread = info::ecf->output.mode == OutputConfiguration::MODE::PTHREAD;

	eslog::info(" ================================ SYSTEM AND ENVIRONMENT INFO ================= %12.3f s\n", eslog::duration());
	eslog::info(" ============================================================================================= \n");

	eslog::info(" == CPU MODEL NAME %72s == \n", cpuinfo.modelName);
	eslog::info(" == NODES %81d == \n", nodes);
	eslog::info(" == SOCKETS PER NODE %70d == \n", cpuinfo.sockets);
	eslog::info(" == CORES PER SOCKET %70d == \n", cpuinfo.cores);
	eslog::info(" == HYPERTHREADING %72s == \n", cpuinfo.hyperthreading ? "ON" : "OFF");
	eslog::info(" == HWTHREADS PER NODE %68d == \n", threads);
	eslog::info(" == SIMD INSTRUCTION SET %66s == \n", simd);
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");

	eslog::info(" == AVAILABLE MEMORY PER NODE [GB] %56.2f == \n", memAvail / 1024 / 1024);
	eslog::info(" == TOTAL MEMORY PER NODE [GB] %60.2f == \n", memTotal / 1024 / 1024);
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");

	eslog::info(" == MPI_COMM_WORLD %72d == \n", info::mpi::size);
	eslog::info(" == OMP_NUM_THREADS %71d == \n", info::env::OMP_NUM_THREADS);
	eslog::info(" == NUMBER OF ASYNCHRONOUS P-THREADS %*d == \n", 54, asynchronousThread);

	eslog::info(" == USED HWTHREADS PER NODE   %*d / %d == \n", 56 - threads / 100, ppn * (info::env::OMP_NUM_THREADS + asynchronousThread), threads);
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");

	eslog::info(" == MPI PROCESSES PER NODE %64d == \n", ppn);
	eslog::info(" == MPI PINNED DOMAIN SIZE %64d == \n", info::system::pinnedDomainSize());
	switch (secondRankNode) {
	case 0: eslog::info(" == MPI PLACEMENT ORDER    %64s == \n", "NODE FIRST"); break;
	case 1: eslog::info(" == MPI PLACEMENT ORDER    %64s == \n", "ACROSS NODES"); break;
	}


	eslog::info(" ============================================================================================= \n\n");
}

void setSignals()
{
	switch (buildmode()) {
	case BUILD_MODE::DEBUG:
	case BUILD_MODE::DEVEL:
		std::signal(SIGTERM, signalHandler);
		std::signal(SIGINT, signalHandler);
	case BUILD_MODE::RELEASE:
	case BUILD_MODE::PROFILE:
		std::signal(SIGSEGV, signalHandler);
		std::signal(SIGILL, signalHandler);
		std::signal(SIGABRT, signalHandler);
		std::signal(SIGFPE, signalHandler);
		break;
	}
}

}
}
}

