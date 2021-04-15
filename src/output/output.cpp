
#include "output.h"

#include "monitors/monitoring.h"
#include "visualization/visualization.h"
#include "visualization/vtklegacy.h"
#include "visualization/ensightgold.h"
#include "visualization/xdmf.h"
#include "visualization/stl.h"
#include "visualization/netgen.h"
#include "visualization/insitu.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/stepinfo.h"
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
	virtual void execute(int tag)
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
};

class AsyncOutputExecutor: public DirectOutputExecutor, public Pthread::Executor, public Pthread {
public:
	AsyncOutputExecutor(): Pthread(this)
	{

	}

	void copy(int tag)
	{
		info::mesh->toBuffer();
		step::toOut();
	}

	void call(int tag)
	{
		DirectOutputExecutor::execute(tag);
	}

	void execute(int tag)
	{
		Pthread::call(tag);
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
		case OutputConfiguration::FORMAT::VTK_LEGACY: writer = new VTKLegacy(info::ecf->output.store_decomposition); break;
		case OutputConfiguration::FORMAT::ENSIGHT: writer = new EnSightGold(info::ecf->output.store_decomposition); break;
		case OutputConfiguration::FORMAT::XDMF: writer = new XDMF(); break;
		case OutputConfiguration::FORMAT::STL_SURFACE: writer = new STL(); break;
		case OutputConfiguration::FORMAT::NETGEN: writer = new Netgen(); break;
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
	if (_async) { _async->execute(OutputExecutor::mesh); }
	if (_direct) { _direct->execute(OutputExecutor::mesh); }
}

void Output::updateMonitors()
{
	if (_allowed) {
		if (_async) { _async->execute(OutputExecutor::monitors); }
		if (_direct) { _direct->execute(OutputExecutor::monitors); }
	}
}

void Output::updateSolution()
{
	if (_allowed) {
		if (_async) { _async->execute(OutputExecutor::solution); }
		if (_direct) { _direct->execute(OutputExecutor::solution); }
	}
}

void Output::suppress()
{
	_allowed = false;
}

void Output::permit()
{
	_allowed = true;
}
