
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "config/ecf/ecf.h"

using namespace espreso;

int main(int argc, char **argv)
{
	info::system::setSignals();
	info::mpi::init(&argc, &argv);

	eslog::init(new Logger<ProgressTerminalLogger>());

	ECF(&argc, &argv);
	eslog::info("ECF syntax is correct\n");

	info::mpi::finish();
	return 0;
}


