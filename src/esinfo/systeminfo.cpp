
#include "systeminfo.h"
#include "eslog.h"
#include "mpi.h"

#include <csignal>
#include <cstdlib>

static void signalHandler(int signal)
{
	switch (signal) {
	case SIGTERM:
		espreso::eslog::error("SIGTERM - termination request, sent to the program.\n");
		break;
	case SIGSEGV:
		espreso::eslog::error("SIGSEGV - invalid memory access (segmentation fault).\n");
		break;
	case SIGINT:
		espreso::eslog::error("SIGINT - external interrupt, usually initiated by the user.\n");
		break;
	case SIGILL:
		espreso::eslog::error("SIGILL - invalid program image, such as invalid instruction.\n");
		break;
	case SIGABRT:
		espreso::eslog::error("SIGABRT - abnormal termination condition, as is e.g. initiated by std::abort().\n");
		break;
	case SIGFPE:
		espreso::eslog::error("SIGFPE - erroneous arithmetic operation such as divide by zero.\n");
		break;
	default:
		espreso::eslog::error("ESPRESO trigger an error.\n");
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

INSTRUCTIONSET instructionSet()
{
	return INSTRUCTIONSET::SSE;
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

