
#ifndef SRC_BASIS_LOGGING_PAPICOUNTERS_H_
#define SRC_BASIS_LOGGING_PAPICOUNTERS_H_

#include "verbosity.h"
#include "esinfo/ecfinfo.h"
#include "wrappers/papi/w.papi.h"

#include <cstddef>
#include <vector>

namespace espreso {

class PAPICounters: public Verbosity<PAPICounters, 'p'> {
public:
	struct Event {
		const char* name;
		long value;

		enum Type {
			START, CHECKPOINT, ACCUMULATED, END, LOADSTEP
		} type;
	};

	struct EventStatistics: public Event {
		long min, max, avg;
		long duration, dmin, dmax, davg;
		long sum, smin, smax, savg; // for accumulated events

		EventStatistics(const Event &event)
		: Event(event),
		  min(event.value), max(event.value), avg(event.value),
		  duration(event.value), dmin(event.value), dmax(event.value), davg(event.value),
		  sum(event.value), smin(event.value), smax(event.value), savg(event.value)
		{ }
	};

	PAPICounters()
	{
		init = 0;
	}

	void initOutput()
	{
		if (info::ecf->output.papi_event.size()) {
			_events.reserve(1000000);
			papi.init();
			init = papi.read();
		}
		if (info::ecf->output.papi_code) {
			_events.reserve(1000000);
			papi.init();
			init = papi.read();
		}
		for (size_t i = 0; i < _events.size(); ++i) {
			_events[i].value = init;
		}
	}

	void start(const char* region, const char* section)
	{
		_events.push_back(Event{ section, papi.read(), Event::START });
	}

	void checkpoint(const char* region)
	{
		_events.push_back(Event{ region, papi.read(), Event::CHECKPOINT });
	}

	void accumulated(const char* region)
	{
		_events.push_back(Event{ region, papi.read(), Event::ACCUMULATED });
	}

	void end(const char* region)
	{
		_events.push_back(Event{ region, papi.read(), Event::END });
	}

	void param(const char* name, const int &value)
	{
		// do nothing
	}

	void param(const char* name, const long &value)
	{
		// do nothing
	}

	void param(const char* name, const long unsigned int &value)
	{
		// do nothing
	}

	void param(const char* name, const double &value)
	{
		// do nothing
	}

	void param(const char* name, const char* value)
	{
		// do nothing
	}

	void ln()
	{
		// do nothing
	}

	void nextLoadStep(int step)
	{
		_events.push_back(Event{ nullptr, papi.read(), Event::LOADSTEP });
	}

	void output(const char* msg, VerboseArg::COLOR color)
	{
		// do nothing
	}

	void error(const char* msg)
	{

	}

	void finish();

protected:
	long init;
	std::vector<Event> _events;
	PAPI papi;
};

}




#endif /* SRC_BASIS_LOGGING_PAPICOUNTERS_H_ */
