
#ifndef SRC_BASIS_LOGGING_TIMELOGGER_H_
#define SRC_BASIS_LOGGING_TIMELOGGER_H_

#include "verbosity.h"

#include <cstddef>
#include <vector>
#include <ctime>

namespace espreso {

class TimeLogger: public Verbosity<TimeLogger, 't'> {
public:
	static double time();
	static double duration();

	struct Event {
		const char* name;

		union Data {
			double time;
			int    ivalue;
			long   lvalue;
			size_t svalue;
			double dvalue;
		} data;

		enum {
			START, CHECKPOINT, ACCUMULATED, END, LOADSTEP,
			INT, LONG, SIZE, DOUBLE
		} type;
	};

	struct EventStatistics: public Event {
		Data min, max, avg;
		Data duration, dmin, dmax, davg;
		Data sum, smin, smax, savg; // for accumulated events

		EventStatistics(const Event &event)
		: Event(event),
		  min(event.data), max(event.data), avg(event.data),
		  duration(event.data), dmin(event.data), dmax(event.data), davg(event.data),
		  sum(event.data), smin(event.data), smax(event.data), savg(event.data)
		{ }
	};

	TimeLogger()
	{
		_events.reserve(1000000);
	}

	void initOutput()
	{

	}

	void start(const char* region, const char* section)
	{
		_events.push_back(Event{ section, Event::Data{ .time = time() }, Event::START });
	}

	void checkpoint(const char* region)
	{
		_events.push_back(Event{ region, Event::Data{ .time = time() }, Event::CHECKPOINT });
	}

	void accumulated(const char* region)
	{
		_events.push_back(Event{ region, Event::Data{ .time = time() }, Event::ACCUMULATED });
	}

	void end(const char* region)
	{
		_events.push_back(Event{ region, Event::Data{ .time = time() }, Event::END });
	}

	void param(const char* name, const int &value)
	{
		_events.push_back(Event{ name, Event::Data{ .ivalue = value }, Event::INT });
	}

	void param(const char* name, const long &value)
	{
		_events.push_back(Event{ name, Event::Data{ .lvalue = value }, Event::LONG });
	}

	void param(const char* name, const long unsigned int &value)
	{
		_events.push_back(Event{ name, Event::Data{ .svalue = value }, Event::SIZE });
	}

	void param(const char* name, const double &value)
	{
		_events.push_back(Event{ name, Event::Data{ .dvalue = value }, Event::DOUBLE });
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
		_events.push_back(Event{ "NEXT STEP", Event::Data{ .ivalue = step }, Event::LOADSTEP });
	}

	void output(const char* msg, VerboseArg::COLOR color)
	{

	}

	void error(const char* msg)
	{

	}

	void finish();

	static time_t initTime;
protected:
	static double initClockTime;
	std::vector<Event> _events;
};

}



#endif /* SRC_BASIS_LOGGING_TIMELOGGER_H_ */
