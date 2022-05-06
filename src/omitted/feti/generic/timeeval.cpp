
#include "timeeval.h"

#include <iomanip>
#include <cmath>
#include <sstream>

#include "mpi.h"
#include "omp.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/eslog.hpp"

using namespace espreso;


// #define SOLVER_MEASUREMENT

// #define MPIBARRIER  ; //MPI_Barrier(info::mpi::MPICommunicator);
#define MPIBARRIER  MPI_Barrier(info::mpi::comm);

#ifndef SOLVER_MEASUREMENT

TimeEvent::TimeEvent(std::string name){};

void TimeEvent::start(){};
void TimeEvent::startWithBarrier(){};
void TimeEvent::startWithoutBarrier(){};

void TimeEvent::start(double time){};
void TimeEvent::startWithBarrier(double time){};
void TimeEvent::startWithoutBarrier(double time){};

void TimeEvent::end(){};
void TimeEvent::endWithBarrier(){};
void TimeEvent::endWithoutBarrier(){};

void TimeEvent::end(double time){};
void TimeEvent::endWithBarrier(double time){};
void TimeEvent::endWithoutBarrier(double time){};

void TimeEvent::reset(){};

void TimeEvent::printStat(double totalTime){};
void TimeEvent::printLastStat(double totalTime){};
double TimeEvent::getLastStat(double totalTime){return 0.0;};

void TimeEvent::printStatMPI(double totalTime){};
void TimeEvent::printLastStatMPI(double totalTime){};
void TimeEvent::printLastStatMPIPerNode(double totalTime){};

void TimeEvent::evaluate(){};
void TimeEvent::evaluateMPI(){};

TimeEval::TimeEval(std::string name):
	totalTime(TimeEvent(name + std::string("- Total "))),
	remainingTime(TimeEvent(name + std::string("- Remaining "))),
	evalName(name)
{
}

void TimeEval::addEvent(TimeEvent &timeEvent){};
void TimeEval::printStats(){};
void TimeEval::printStatsMPI(){};

#else

TimeEvent::TimeEvent(std::string name)
{
	eventName   = name;
	name_length = 50;
	val_length  = 12;
	reset();
}



void TimeEvent::reset() {
	eventCount = 0;
	eventTime.clear();

	avgTime = 0.0;
	sumTime = 0.0;
	minTime = 0.0;
	maxTime = 0.0;
	stdDev  = 0.0;

	g_avgTime = 0.0;
	g_minTime = 0.0;
	g_maxTime = 0.0;
	g_stdDev  = 0.0;
	g_sumTime = 0.0;
}

double TimeEvent::time()
{
	return omp_get_wtime();
}

void TimeEvent::start(double time) {
#ifdef TM_BLOCK_START
	startWithBarrier(time);
#else
	startWithoutBarrier(time);
#endif
}

void TimeEvent::startWithBarrier(double time) {
	MPIBARRIER;
	startWithoutBarrier(time);
}

void TimeEvent::startWithoutBarrier(double time) {
	eventTime.push_back(time);
}

void TimeEvent::start() {
#ifdef TM_BLOCK_START
	startWithBarrier();
#else
	startWithoutBarrier();
#endif
}

void TimeEvent::startWithBarrier() {
	MPIBARRIER;
	startWithoutBarrier();
}

void TimeEvent::startWithoutBarrier() {
	eventTime.push_back(time());
}


void TimeEvent::end() {
#ifdef TM_BLOCK_END
	endWithBarrier();
#else
	endWithoutBarrier();
#endif
}

void TimeEvent::endWithBarrier() {
	MPIBARRIER;
	endWithoutBarrier();
}

void TimeEvent::endWithoutBarrier() {
	eventTime.back() = time() - eventTime.back();
	eventCount++;
}

void TimeEvent::end(double time) {
#ifdef TM_BLOCK_END
	endWithBarrier(time);
#else
	endWithoutBarrier(time);
#endif
}

void TimeEvent::endWithBarrier(double time) {
	MPIBARRIER;
	endWithoutBarrier(time);
}

void TimeEvent::endWithoutBarrier(double time) {
	eventTime.back() = time - eventTime.back();
	eventCount++;
}


void TimeEvent::evaluate() {
	sumTime = 0;
	avgTime = 0;
	minTime = eventCount ? eventTime[0] : 0;
	maxTime = 0;
	stdDev  = 0;

	for (esint i = 0; i < eventCount; i++) {
		double d_time = eventTime[i];
		sumTime += d_time;
		if (d_time < minTime) {
			minTime = d_time;
		}
		if (d_time > maxTime) {
			maxTime = d_time;
		}
	}

	avgTime = sumTime / eventCount;

	double E = 0;
	for (esint i = 0; i < eventCount; i++) {
		E += (eventTime[i] - avgTime) * (eventTime[i] - avgTime);
	}

	if (eventCount * E > 0) {
		stdDev = sqrt(1 / eventCount * E);
	} else {
		stdDev = 0;
	}
}

void TimeEvent::evaluateMPI() {
	evaluate();

	MPI_Reduce(&avgTime, &g_avgTime, 1, MPI_DOUBLE, MPI_SUM, 0, info::mpi::comm);
	g_avgTime = g_avgTime / info::mpi::size;

	MPI_Reduce(&sumTime, &g_sumTime, 1, MPI_DOUBLE, MPI_SUM, 0, info::mpi::comm);
	g_sumTime = g_sumTime / info::mpi::size;

	MPI_Reduce(&minTime, &g_minTime, 1, MPI_DOUBLE, MPI_MIN, 0, info::mpi::comm);
	MPI_Reduce(&maxTime, &g_maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, info::mpi::comm);
}


void TimeEvent::printStat(double totalTime) {
	evaluate();

	std::stringstream ss;
	ss	<< std::setw(name_length) << std::left << eventName
		<< " avg.: " << std::fixed << std::setw(val_length) << avgTime
		<< " min.: " << std::setw(val_length) << minTime
		<< " max.: " << std::setw(val_length) << maxTime
		<< " % of avg tot: " << std::setw(val_length)
		<< (totalTime != 0 ? 100.0 * avgTime / totalTime : INFINITY);
	eslog::duration("%s\n", ss.str().c_str());
}


void TimeEvent::printLastStat(double totalTime) {
	avgTime = eventTime.back();

	std::stringstream ss;
	ss	<< std::setw(name_length) << std::left << eventName
		<< " avg.: " << std::fixed << std::setw(val_length) << avgTime
		<< " min.: " << std::setw(val_length) << " -- "
		<< " max.: " << std::setw(val_length) << " -- "
		<< " % of avg tot: " << std::setw(val_length)
		<< (totalTime != 0 ? 100.0 * avgTime / totalTime : INFINITY);
	eslog::duration("%s\n", ss.str().c_str());
}

double TimeEvent::getLastStat(double totalTime) {
	return eventTime.back();
}


void TimeEvent::printStatMPI(double totalTime) {
	evaluateMPI();
	std::stringstream ss;
	ss	<< std::setw(name_length) << std::left << eventName
		<< " avg.: " << std::setw(val_length) << std::fixed << g_avgTime
		<< " min.: " << std::setw(val_length) << g_minTime
		<< " max.: " << std::setw(val_length) << g_maxTime
		<< " sum.: " << std::setw(val_length) << g_sumTime
		<< " count: " << std::setw(val_length) << eventCount
		<< " % of avg tot: " << std::setw(val_length)
		<< (totalTime != 0 ? 100.0 * g_avgTime / totalTime : INFINITY);
	eslog::duration("%s\n", ss.str().c_str());
}


void TimeEvent::printLastStatMPI(double totalTime) {
	if (info::ecf->output.measure_level == 0) {
		return;
	}
	double d_time = eventTime.back();

	MPI_Reduce(&d_time, &g_avgTime, 1, MPI_DOUBLE, MPI_SUM, 0, info::mpi::comm);
	g_avgTime = g_avgTime / info::mpi::size;

	MPI_Reduce(&d_time, &g_minTime, 1, MPI_DOUBLE, MPI_MIN, 0, info::mpi::comm);
	MPI_Reduce(&d_time, &g_maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, info::mpi::comm);

	std::stringstream ss;
	ss	<< std::setw(name_length) << std::left << eventName
		<< " last: " << std::fixed << std::setw(val_length) << g_avgTime
		<< " min.: " << std::setw(val_length) << g_minTime
		<< " max.: " << std::setw(val_length) << g_maxTime
		<< " % of avg tot: " << std::setw(val_length)
		<< (totalTime != 0 ? 100.0 * g_avgTime / totalTime : INFINITY);
	eslog::duration("%s\n", ss.str().c_str());
}


void TimeEvent::printLastStatMPIPerNode(double totalTime)
{
	double d_time = eventTime.back();
	std::vector<double> d_all_times;

	if(info::mpi::rank == 0) {
		d_all_times.resize(info::mpi::size);
	} else {
		d_all_times.resize(1);
	}

	MPI_Gather(&d_time, 1, MPI_DOUBLE, &d_all_times[0], 1, MPI_DOUBLE, 0, info::mpi::comm);

	MPI_Reduce(&d_time, &g_avgTime, 1, MPI_DOUBLE, MPI_SUM, 0, info::mpi::comm);
	g_avgTime= g_avgTime / info::mpi::size;

	MPI_Reduce(&d_time, &g_minTime, 1, MPI_DOUBLE, MPI_MIN, 0, info::mpi::comm);
	MPI_Reduce(&d_time, &g_maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, info::mpi::comm);

	std::stringstream ss;
	ss	<< std::setw(name_length) << std::left << eventName
		<< " last: " << std::fixed << std::setw(val_length) << g_avgTime
		<< " min.: " << std::setw(val_length) << g_minTime
		<< " max.: " << std::setw(val_length) << g_maxTime
		<< " % of avg tot: " << std::setw(val_length)
		<< (totalTime != 0 ? 100.0 * g_avgTime / totalTime : INFINITY);
	eslog::duration("%s\n", ss.str().c_str());


	ss.clear();
	for (esint i = 0; i < info::mpi::size; i++) {
		ss << std::fixed << std::setw(3) << "R: " << std::setw(5) << i << std::setw(15) << d_all_times[i];

		if ((i + 1) % 10 == 0) {
			eslog::duration("%s\n", ss.str().c_str());
			ss.clear();
		}
	}
	eslog::duration("%s\n", ss.str().c_str());
}



TimeEval::TimeEval(std::string name):
	totalTime(TimeEvent(name + std::string("- Total "))),
	remainingTime(TimeEvent(name + std::string("- Remaining "))),
	evalName(name)
{
}

void TimeEval::addEvent(TimeEvent &timeEvent) {
	timeEvents.push_back(timeEvent);
}

void TimeEval::addPointerToEvent(TimeEvent *timeEvent)
{
	ptimeEvents.push_back(timeEvent);
}

void TimeEval::printStats() {
	totalTime.evaluate();

	for (size_t i = 0; i < timeEvents.size(); i++) {
		timeEvents[i].printStat(totalTime.avgTime);
	}

	totalTime.printStat(totalTime.avgTime);
}

void TimeEval::printStatsMPI() {
	double sum_avg_time = 0;

	auto separator = [] (int size, char character) {
		std::stringstream ss;
		for (int i = 0; i < size; i++) {
			ss << character;
		}
		return ss.str();
	};

	int separator_size = 80;

	std::stringstream ss;
	ss << separator(separator_size, '*') << "\n";
	ss << "        " << evalName << "         \n";
	ss << separator(separator_size, '*') << "\n";
	eslog::duration("%s", ss.str().c_str());
	totalTime.evaluateMPI();

	for (size_t i = 0; i < timeEvents.size(); i++) {
		timeEvents[i].printStatMPI(totalTime.g_avgTime);
		sum_avg_time += timeEvents[i].avgTime;
	}

	for (size_t i = 0; i < ptimeEvents.size(); i++) {
		ptimeEvents[i]->printStatMPI(totalTime.g_avgTime);
		sum_avg_time += ptimeEvents[i]->avgTime;
	}

	{
		std::stringstream sss;
		sss << separator(separator_size, '-');
		eslog::duration("%s\n", sss.str().c_str());
	}

	totalTime.printStatMPI(totalTime.g_avgTime);

	{
		std::stringstream sss;
		sss << separator(separator_size, '*');
		eslog::duration("%s\n", ss.str().c_str());
	}
}

#endif




