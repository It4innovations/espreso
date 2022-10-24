
#include "output.h"

#include "monitors/monitoring.h"
#include "visualization/visualization.h"
#include "visualization/vtklegacy.h"
#include "visualization/ensightgold.h"
#include "visualization/ensightgold.volume.h"
#include "visualization/xdmf.h"
#include "visualization/stl.h"
#include "visualization/netgen.h"
#include "visualization/insitu.h"
#include "visualization/openvdb.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/stepinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "wrappers/pthread/w.pthread.h"

#include <vector>

namespace espreso {

struct OutputExecutor {
	static const int mesh = 0, solution = 1, monitors = 2;

	void insert(OutputWriter *writer)
	{
		writers.push_back(writer);
	}

	virtual ~OutputExecutor()
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			delete writers[i];
		}
	}

	virtual void execute(int tag) = 0;

	std::vector<OutputWriter*> writers;
};

class DirectOutputExecutor: public OutputExecutor {
public:
	void call(int tag)
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			if (tag == mesh) {
				writers[i]->updateMesh();
			}
			if (tag == solution) {
				writers[i]->updateSolution();
			}
			if (tag == monitors) {
				writers[i]->updateMonitors();
			}
		}
	}

	virtual void execute(int tag)
	{
		step::toOut();
		if (tag == mesh || tag == solution) {
			for (size_t i = 0; i < writers.size(); ++i) {
				writers[i]->lock();
			}
		}
		call(tag);
	}
};

class AsyncOutputExecutor: public DirectOutputExecutor, public Async {
public:
	AsyncOutputExecutor(): Async(this)
	{

	}

	void copy(int tag)
	{
		info::mesh->toBuffer();
		step::toOut();
		if (tag == mesh || tag == solution) {
			for (size_t i = 0; i < writers.size(); ++i) {
				writers[i]->lock();
			}
		}
	}

	void call(int tag)
	{
		DirectOutputExecutor::call(tag);
	}

	void execute(int tag)
	{
		Pthread::execute(tag);
	}
};

}

using namespace espreso;

OutputWriter::OutputWriter()
: _path(info::ecf->outpath + "/"), _directory("PREPOSTDATA/"), _name(info::ecf->name),
  _measure(info::ecf->output.mode == OutputConfiguration::MODE::SYNC), _allowed(true)
{

}

void OutputWriter::createOutputDirectory()
{
	utils::createDirectory({ _path, _directory });
}

Output::Output()
: _direct(new DirectOutputExecutor()), _async(NULL)
{
	if (info::ecf->output.mode != OutputConfiguration::MODE::SYNC) {
		_async = new AsyncOutputExecutor();
	}
	if (info::ecf->output.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		OutputWriter *writer = NULL;
		switch (info::ecf->output.format) {
		case OutputConfiguration::FORMAT::VTK_LEGACY: writer = new VTKLegacy(); break;
		case OutputConfiguration::FORMAT::ENSIGHT: writer = new EnSightGold(); break;
		case OutputConfiguration::FORMAT::ENSIGHT_VOLUME: writer = new EnSightGoldVolume(); break;
		case OutputConfiguration::FORMAT::XDMF: writer = new XDMF(); break;
		case OutputConfiguration::FORMAT::STL_SURFACE: writer = new STL(); break;
		case OutputConfiguration::FORMAT::NETGEN: writer = new Netgen(); break;
		case OutputConfiguration::FORMAT::OPENVDB: writer = new OpenVDB(); break;
		default:
			eslog::internalFailure("implement the selected output format.\n");
		}
		switch (info::ecf->output.mode) {
		case OutputConfiguration::MODE::SYNC: _direct->insert(writer); break;
		case OutputConfiguration::MODE::PTHREAD: _async->insert(writer); break;
		default:
			eslog::internalFailure("implement the selected output mode.\n");
		}
	}
	if (info::ecf->output.monitors_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		_direct->insert(new Monitoring());
	}
}

Output::~Output()
{
	if (_direct) { delete _direct; }
	if (_async) { delete _async; }
}

void Output::updateMesh()
{
	eslog::info(" ===================================== STORE GEOMETRY ========================= %12.3f s\n", eslog::duration());
	switch (info::ecf->output.format) {
	case OutputConfiguration::FORMAT::ENSIGHT:        eslog::info(" == OUTPUT FORMAT %73s == \n", "ENSIGHT"); break;
	case OutputConfiguration::FORMAT::ENSIGHT_VOLUME: eslog::info(" == OUTPUT FORMAT %73s == \n", "ENSIGHT VOLUME"); break;
	case OutputConfiguration::FORMAT::NETGEN:         eslog::info(" == OUTPUT FORMAT %73s == \n", "NETGET"); break;
	case OutputConfiguration::FORMAT::OPENVDB:        eslog::info(" == OUTPUT FORMAT %73s == \n", "OPEN VDB"); break;
	case OutputConfiguration::FORMAT::STL_SURFACE:    eslog::info(" == OUTPUT FORMAT %73s == \n", "STL SURFACE"); break;
	case OutputConfiguration::FORMAT::VTK_LEGACY:     eslog::info(" == OUTPUT FORMAT %73s == \n", "VTK LEGACY"); break;
	case OutputConfiguration::FORMAT::XDMF:           eslog::info(" == OUTPUT FORMAT %73s == \n", "XDMF"); break;
	}

	switch (info::ecf->output.mode) {
	case OutputConfiguration::MODE::SYNC: eslog::info(" == OUTPUT MODE %75s == \n", "SYNC"); break;
	case OutputConfiguration::MODE::PTHREAD: eslog::info(" == OUTPUT MODE %75s == \n", "PTHREAD"); break;
	}
	if (_async) { _async->execute(OutputExecutor::mesh); }
	if (_direct) { _direct->execute(OutputExecutor::mesh); }
	eslog::info(" ============================================================================================= \n\n");
	eslog::checkpointln("ESPRESO: MESH STORED");
}

void Output::updateMonitors()
{
	if (_allowed) {
		if (_async) { _async->execute(OutputExecutor::monitors); }
		if (_direct) { _direct->execute(OutputExecutor::monitors); }
	}
	eslog::checkpointln("MESH: MONITORS UPDATED");
}

void Output::updateSolution()
{
	if (_allowed) {
		if (_async) { _async->execute(OutputExecutor::solution); }
		if (_direct) { _direct->execute(OutputExecutor::solution); }
	}
	eslog::checkpointln("MESH: SOLUTION UPDATED");
}

void Output::suppress()
{
	_allowed = false;
}

void Output::permit()
{
	_allowed = true;
}
