
#ifndef BASIS_LOGGING_LOGGING_H_
#define BASIS_LOGGING_LOGGING_H_

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sys/sysinfo.h>

#include "esconfig.h"
#include "timeeval.h"

#define ESTEST(EVENT) if (!espreso::Test::report(EVENT))    ; else espreso::Test(EVENT).get()
#define ESINFO(EVENT) if (!espreso::Info::report(EVENT))    ; else espreso::Info(EVENT).get()
#define ESLOG(EVENT)  if (!espreso::Measure::report(EVENT)) ; else espreso::Measure(EVENT).get()

namespace espreso {

enum ESPRESOTest {
	FAILED,
	PASSED
};

enum TestEvent {
	MANDATORY,
	TEST_LEVEL0,

	SIMPLE,
	TEST_LEVEL1,

	EXPENSIVE,
	TEST_LEVEL2,

	PEDANTIC,
	TEST_LEVEL3,
};

enum InfoEvent {
	ERROR,
	ALWAYS,
	VERBOSE_LEVEL0,

	CONVERGENCE,
	OVERVIEW,
	PROGRESS1,
	VERBOSE_LEVEL1,

	DETAILS,
	PROGRESS2,
	VERBOSE_LEVEL2,

	EXHAUSTIVE,
	PROGRESS3,
	LIBRARIES,
	VERBOSE_LEVEL3
};

enum MeasureEvent {
	MEASURE_LEVEL0,

	SUMMARY,
	CHECKPOINT1,
	MEASURE_LEVEL1,

	CHECKPOINT2,
	MEASURE_LEVEL2,

	MEMORY,
	CHECKPOINT3,
	MEASURE_LEVEL3
};

class Test
{
public:
	Test& operator<<(const ESPRESOTest &test)
	{
		if (test == FAILED) { error = true; }
		return *this;
	}
	template<typename Ttype>
	Test& operator<<(const Ttype &value)
	{
		os << value;
		return *this;
	}

	Test(TestEvent event): error(false) {};
	~Test();

	Test& get() { return *this; };

	static bool report(TestEvent event) {
		switch (config::info::testingLevel) {
		case 0: return event < TEST_LEVEL0;
		case 1: return event < TEST_LEVEL1;
		case 2: return event < TEST_LEVEL2;
		case 3: return event < TEST_LEVEL3;
		default : return true;
		}
	};

protected:
	std::ostringstream os;
	bool error;
};

class Info
{
public:
	Info(InfoEvent event): event(event), _plain(false) {};
	~Info();

	enum InfoMode { FORMATTED, PLAIN };
	Info& operator<<(const InfoMode& mode)
	{
		_plain = mode == PLAIN;
		return *this;
	}
	template<typename Ttype>
	Info& operator<<(const Ttype &value)
	{
		os << value;
		return *this;
	}

	Info& get() { return *this; };

	static InfoMode plain() { return PLAIN; }

	template<typename Tvalue>
	static std::string sumValue(const Tvalue &value);
	template<typename Tvalue>
	static std::string averageValue(const Tvalue &value);


	static bool report(InfoEvent event) {
		switch (config::info::verboseLevel) {
		case 0: return event < VERBOSE_LEVEL0;
		case 1: return event < VERBOSE_LEVEL1;
		case 2: return event < VERBOSE_LEVEL2;
		case 3: return event < VERBOSE_LEVEL3;
		default : return true;
		}
	};

protected:
	std::ostringstream os;
	InfoEvent event;
	bool _plain;
};

struct Checkpoint;

class Measure
{
public:
	static double time()
	{
		return omp_get_wtime();
	}

	Measure(MeasureEvent event): event(event) {};
	~Measure();

	static double processMemory();
	static double usedRAM();
	static double availableRAM();

	std::ostringstream& get() { return os; };

	static bool report(MeasureEvent event) {
		switch (config::info::measureLevel) {
		case 0: return event < MEASURE_LEVEL0;
		case 1: return event < MEASURE_LEVEL1;
		case 2: return event < MEASURE_LEVEL2;
		case 3: return event < MEASURE_LEVEL3;
		default : return true;
		}
	};

protected:
	static std::vector<Checkpoint> checkpoints;

	void evaluateCheckpoints();

	std::ostringstream os;
	MeasureEvent event;
};

class Logging {

public:
	static std::string prepareFile(const std::string &name)
	{
		std::stringstream dir, file, mkdir;

		dir << config::info::output << "/" << config::MPIrank << "/";
		file << dir.str() << "/" << name << ".txt";

		mkdir << "mkdir -p " << dir.str();
		system(mkdir.str().c_str());

		return file.str();
	}

	static std::string prepareFile(size_t subdomain, const std::string &name)
	{
		std::stringstream ss;
		ss << name << subdomain;
		return prepareFile(ss.str());
	}
};
}

#include "logging.hpp"


#endif /* BASIS_LOGGING_LOGGING_H_ */
