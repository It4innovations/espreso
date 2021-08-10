
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

	virtual void mesh() = 0;
	virtual void monitors(step::TYPE type) = 0;
	virtual void solution(const step::Time &time) = 0;
	virtual void solution(const step::Frequency &frequency) = 0;

	std::vector<OutputWriter*> writers;
};

class DirectOutputExecutor: public OutputExecutor {
public:
	virtual void mesh()
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			writers[i]->updateMesh();
		}
	}

	virtual void monitors(step::TYPE type)
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			writers[i]->updateMonitors(type);
		}
	}

	virtual void solution(const step::Time &time)
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			writers[i]->updateSolution(time);
		}
	}

	virtual void solution(const step::Frequency &frequency)
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			writers[i]->updateSolution(frequency);
		}
	}
};



class AsyncOutputExecutor: public DirectOutputExecutor, public Pthread::Executor, public Pthread {
	struct SharedData {
		enum class TAG {
			MESH,
			MONITORS,
			SOLUTION,
			DUMMY
		} tag;
		step::TYPE type;
		step::Time time;
		step::Frequency frequency;

		SharedData(TAG tag, step::TYPE type): tag(tag), type(type) {}
		SharedData(TAG tag, const step::Time &time): tag(tag), type(step::TYPE::TIME), time(time) {}
		SharedData(TAG tag, const step::Frequency &frequency): tag(tag), type(step::TYPE::FREQUENCY), frequency(frequency) {}
	} app, thread;

public:
	AsyncOutputExecutor(): Pthread(this), app(SharedData::TAG::DUMMY, step::TYPE::TIME), thread(app)
	{

	}

	void copy()
	{
		info::mesh->toBuffer();
		thread = app;
	}

	void call()
	{
		if (thread.tag == SharedData::TAG::MESH) {
			DirectOutputExecutor::mesh();
		}
		if (thread.tag == SharedData::TAG::MONITORS) {
			DirectOutputExecutor::monitors(thread.type);
		}
		if (thread.tag == SharedData::TAG::SOLUTION && thread.type == step::TYPE::TIME) {
			DirectOutputExecutor::solution(thread.time);
		}
		if (thread.tag == SharedData::TAG::SOLUTION && thread.type == step::TYPE::FREQUENCY) {
			DirectOutputExecutor::solution(thread.frequency);
		}
	}

	virtual void mesh()
	{
		app = SharedData(SharedData::TAG::MESH, step::TYPE::TIME); // dummy type
		Pthread::call();
	}

	virtual void monitors(step::TYPE type)
	{
		app = SharedData(SharedData::TAG::MONITORS, type);
		Pthread::call();
	}

	virtual void solution(const step::Time &time)
	{
		app = SharedData(SharedData::TAG::SOLUTION, time);
		Pthread::call();
	}

	virtual void solution(const step::Frequency &frequency)
	{
		app = SharedData(SharedData::TAG::SOLUTION, frequency);
		Pthread::call();
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
	// Time is not used during storing the solution
	if (_async) { _async->mesh(); }
	if (_direct) { _direct->mesh(); }
}

void Output::updateMonitors(step::TYPE type)
{
	if (_allowed) {
		if (_async) { _async->monitors(type); }
		if (_direct) { _direct->monitors(type); }
	}
}

void Output::updateSolution(const step::Time &time)
{
	if (_allowed) {
		if (_async) { _async->solution(time); }
		if (_direct) { _direct->solution(time); }
	}
}

void Output::updateSolution(const step::Frequency &frequency)
{
	if (_allowed) {
		if (_async) { _async->solution(frequency); }
		if (_direct) { _direct->solution(frequency); }
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
