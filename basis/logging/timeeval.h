
#ifndef BASIS_LOGGING_TIMEEVAL_H_
#define BASIS_LOGGING_TIMEEVAL_H_

#include "mpi.h"
#include "omp.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

struct Checkpoint
{
	Checkpoint(const std::string &name, double time, size_t level)
	:name(name), time(time), level(level) { };

	double time;
	size_t level;
	std::string name;
};

class Interval
{
public:
	Interval(const std::string &name, double duration): name(name), duration(duration) {};

	void insert(const Interval &interval)
	{
		subIntervals.push_back(interval);
	}
protected:
	double duration;
	std::string name;

	std::vector<Interval> subIntervals;
};

struct TimeEvent
{
	friend class TimeEval;

	TimeEvent(std::string name);

	static double time()
	{
		return omp_get_wtime();
	}

	void start();
	void startWithBarrier();
	void startWithoutBarrier();
	void start(double time);
	void startWithBarrier(double time);
	void startWithoutBarrier(double time);

	void end();
	void endWithBarrier();
	void endWithoutBarrier();
	void end(double time);
	void endWithBarrier(double time);
	void endWithoutBarrier(double time);

	void reset();

	void printStat(double totalTime = 0.0);
	void printLastStat(double totalTime = 0.0);

	void printStatMPI(double totalTime = 0.0);
	void printLastStatMPI(double totalTime = 0.0);
	void printLastStatMPIPerNode(double totalTime = 0.0);

private:
	void evaluate();
	void evaluateMPI();

	std::string eventName;
	eslocal eventCount;
	std::vector<double> eventTime;

	eslocal name_length;
	eslocal val_length;

	double avgTime;
	double sumTime;
	double minTime;
	double maxTime;
	double stdDev;

	double g_avgTime;
	double g_sumTime;
	double g_minTime;
	double g_maxTime;
	double g_stdDev;
};


struct TimeEval
{
	TimeEval(std::string name);

	void addEvent(TimeEvent &timeEvent);
	void printStats();
	void printStatsMPI();

	TimeEvent totalTime;
	TimeEvent remainingTime;

	std::string evalName;

	std::string eventName;
	std::vector<TimeEvent> timeEvents;

};

#endif /* BASIS_LOGGING_TIMEEVAL_H_ */
