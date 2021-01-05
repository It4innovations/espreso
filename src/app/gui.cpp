
#include "wrappers/mpi/communication.h"

#include "esinfo/eslog.hpp"
#include "esinfo/stepinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"

#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/timelogger.h"
#include "basis/logging/profiler.h"
#include "basis/utilities/sysutils.h"

#include "config/reader/reader.h"
#include "config/configuration.h"
#include "mesh/mesh.h"
#include "output/resultstore.h"
#include "physics/physicalsolver.h"

#include "gui/parallel/mpimanager.h"
#include "gui/workspace/workspacewindow.h"

#include <locale>
#include <QDebug>
#include <QApplication>

using namespace espreso;

int main(int argc, char *argv[])
{
	profiler::syncstart("espreso");

	profiler::syncstart("initialization");
	info::system::setSignals();
	info::env::set();
	profiler::synccheckpoint("set_signals_and_env");
	info::mpi::init(&argc, &argv);
	profiler::synccheckpoint("mpi_init");
	MPITools::init();
	profiler::synccheckpoint("mpi_init_tools");

	eslog::init(new Logger<TimeLogger, ProgressTerminalLogger, ProgressFileLogger>);
	profiler::synccheckpoint("init_loggers");
	eslog::startln("ESPRESO: STARTED", "ESPRESO");

	ECF::init(&argc, &argv, "espresogui");
	profiler::synccheckpoint("init_configuration");
	eslog::checkpointln("ESPRESO: CONFIGURATION READ");
	eslog::startln("CONFIGURATION STARTED", "CONFIGURATION");

	bool divided = info::mpi::divide(info::ecf->input.decomposition.mesh_duplication);
	MPITools::setSubset(info::mpi::size); //info::ecf->input.third_party_scalability_limit);
	eslog::initFiles();
	profiler::synccheckpoint("divide_mpi");
	eslog::printRunInfo(&argc, &argv);
	profiler::synccheckpoint("init_run_info");
	if (!divided) {
		eslog::globalerror("Cannot set MESH DUPLICATION: the number of MPI processes is not divisible by %d\n", info::ecf->input.decomposition.mesh_duplication);
	}
	eslog::checkpointln("CONFIGURATION: RUN INFO INITIALIZED");

	Mesh::init();
	profiler::synccheckpoint("init_mesh");
	eslog::endln("CONFIGURATION: MESH INITIALIZED");
	eslog::checkpointln("ESPRESO: RUN INITIALIZED");
	profiler::syncend("initialization");

	{ // run GUI
		QApplication gui(argc, argv);

		MpiManager mpim(argc, argv);

		WorkspaceWindow ww(&mpim);

		setlocale(LC_ALL, "C");
		setlocale(LC_CTYPE, "C");

		if (info::mpi::rank == 0)
		{
			ww.init();
			ww.showMaximized();

			// CSS
	//        QCoreApplication::setAttribute(Qt::AA_UseStyleSheetPropagationInWidgetStyles, true);
	//        QFile css(":/stylesheets/stylesheets/general.css");
	//        css.open(QFile::ReadOnly);
	//        QTextStream stream(&css);
	//        a.setStyleSheet(stream.readAll());
		}

		mpim.loop();

		if (info::mpi::rank == 0)
		{
			gui.exec();
			mpim.masterExit();
		}
	}

	eslog::finish();
	profiler::syncend("espreso");
	profiler::print(); // need to be printed before MPI_Finalize

	Mesh::finish();
	ECF::finish();
	MPITools::finish();
	info::mpi::finish();
	return 0;
}
