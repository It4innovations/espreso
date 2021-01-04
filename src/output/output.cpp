
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
#include "mesh/mesh.h"
#include "wrappers/pthread/w.pthread.h"

namespace espreso {

struct OutputExecutor {
	static const int mesh = 0, solution = 1, monitors = 2;

	OutputExecutor(OutputWriter * writer): writer(writer) {}
	virtual ~OutputExecutor() { delete writer; }

	virtual void execute(int tag) = 0;

	OutputWriter *writer;
};

class DirectOutputExecutor: public OutputExecutor {
public:
	DirectOutputExecutor(OutputWriter * writer): OutputExecutor(writer) {}

	virtual void execute(int tag)
	{
		if (tag == mesh) {
			writer->updateMesh();
		}
		if (tag == solution) {
			writer->updateSolution();
		}
		if (tag == monitors) {
			writer->updateMonitors();
		}
	}
};

class AsyncOutputExecutor: public DirectOutputExecutor, public Pthread, public Pthread::Executor {
public:
	AsyncOutputExecutor(OutputWriter * writer): DirectOutputExecutor(writer), Pthread(this)
	{

	}

	void call(int tag)
	{
		DirectOutputExecutor::execute(tag);
	}

	void execute(int tag)
	{
		if (tag == monitors) {
			DirectOutputExecutor::execute(tag);
		} else {
			Pthread::call(tag);
		}
	}
};

}

using namespace espreso;
OutputWriter::OutputWriter()
: _path(info::ecf->outpath + "/"), _directory("PREPOSTDATA/"), _name(info::ecf->name),
  _measure(Mesh::convertDatabase()), _allowed(true)
{

}

void OutputWriter::createOutputDirectory()
{
	utils::createDirectory({ _path, _directory });
}

template <class TExecutor>
static OutputExecutor* createWriter()
{
	switch (info::ecf->output.format) {
	case OutputConfiguration::FORMAT::VTK_LEGACY:
		return new TExecutor(new VTKLegacy(!info::ecf->input.convert_database));
	case OutputConfiguration::FORMAT::ENSIGHT:
		return new TExecutor(new EnSightGold(!info::ecf->input.convert_database));
	case OutputConfiguration::FORMAT::XDMF:
		return new TExecutor(new XDMF());
	case OutputConfiguration::FORMAT::STL_SURFACE:
		return new TExecutor(new STL());
	case OutputConfiguration::FORMAT::NETGEN:
		return new TExecutor(new Netgen());
	default:
		eslog::internalFailure("implement the selected output format.\n");
		return NULL;
	}
}

Output::Output()
{
	if (info::ecf->output.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
		switch (info::ecf->output.mode) {
		case OutputConfiguration::MODE::SYNC: _executors.push_back(createWriter<DirectOutputExecutor>()); break;
		case OutputConfiguration::MODE::PTHREAD: _executors.push_back(createWriter<AsyncOutputExecutor>()); break;
		}
	}
}

Output::~Output()
{
	for (size_t i = 0; i < _executors.size(); ++i) {
		delete _executors[i];
	}
}

void Output::updateMesh()
{
	for (size_t i = 0; _allowed && i < _executors.size(); ++i) {
		_executors[i]->execute(OutputExecutor::mesh);
	}
}

void Output::updateMonitors()
{
	for (size_t i = 0; _allowed && i < _executors.size(); ++i) {
		if (_executors[i]->writer->storeStep()) {
			_executors[i]->execute(OutputExecutor::monitors);
		}
	}
}

void Output::updateSolution()
{
	for (size_t i = 0; _allowed && i < _executors.size(); ++i) {
		if (_executors[i]->writer->storeStep()) {
			_executors[i]->execute(OutputExecutor::solution);
		}
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
