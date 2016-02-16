
#ifndef BASIS_LOGGING_LOGGING_H_
#define BASIS_LOGGING_LOGGING_H_

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>

#include "esconfig.h"
#include "timeeval.h"

#define ESLOG(EVENT) if (!eslog::Log::report(EVENT)) ; else eslog::Log(EVENT).get()
#define ESTEST(EVENT) if (!eslog::Test::report(EVENT)) ; else eslog::Test(EVENT).get()

namespace eslog {

enum ESPRESOTest {
	FAILED,
	PASSED
};

enum TestEvent {
	SIMPLE,
	EXPENSIVE,
	PEDANTIC
};

enum LogEvent {
	ERROR,
	VERBOSE_LEVEL0,
	INFO,
	CHECKPOINT1,
	SUMMARY,
	VERBOSE_LEVEL1,
	CHECKPOINT2,
	VERBOSE_LEVEL2,
	CHECKPOINT3,
	VERBOSE_LEVEL3
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
		return event < esconfig::info::testingLevel;
	};

protected:
	std::ostringstream os;
	bool error;
};

class Log
{
public:
	static double time()
	{
		return omp_get_wtime();
	}

	Log(LogEvent event);
	~Log();

	std::ostringstream& get() { return os; };

	static bool report(LogEvent event) {
		switch (esconfig::info::verboseLevel) {
		case 0: return event < VERBOSE_LEVEL0;
		case 1: return event < VERBOSE_LEVEL1;
		case 2: return event < VERBOSE_LEVEL2;
		case 3: return event < VERBOSE_LEVEL3;
		default : return true;
		}
	};

protected:
	static std::vector<Checkpoint> checkpoints;

	void evaluateCheckpoints();

	std::ostringstream os;
	LogEvent event;
};



class Logging {

public:
	static std::string prepareFile(const std::string &name)
	{
		std::stringstream dir, file, mkdir;

		dir << esconfig::info::output << "/" << esconfig::MPIrank << "/";
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


#endif /* BASIS_LOGGING_LOGGING_H_ */
