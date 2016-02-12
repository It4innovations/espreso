
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

namespace eslog {

enum ESPRESOTest {
	FAILED,
	PASSED
};

enum Event {
	ERROR,
	VERBOSE_LEVEL0,
	DURATION,
	CHECKPOINT1,
	CHECKPOINT2,
	CHECKPOINT3,
	TEST_SIMPLE,
	VERBOSE_LEVEL1,
	VERBOSE_LEVEL2,
	TEST_EXPENSIVE,
	VERBOSE_LEVEL3
};

class Log
{
public:
	Log& operator<<(const ESPRESOTest &test)
	{
		if (test == FAILED) { error = true; }
		return *this;
	}
	template<typename Ttype>
	Log& operator<<(const Ttype &value)
	{
		os << value;
		return *this;
	}

	Log(Event event);
	~Log();

	Log& get() { return *this; };

	static bool report(Event event) {
		switch (esconfig::info::verboseLevel) {
		case 0: return event < VERBOSE_LEVEL0;
		case 1: return event < VERBOSE_LEVEL1;
		case 2: return event < VERBOSE_LEVEL2;
		case 3: return event < VERBOSE_LEVEL3;
		default : return true;
		}
	};

protected:
	static std::vector<double> lastTimes;
	static double start;

	std::ostringstream os;
	Event event;
	bool error;
private:
	Log(const Log&);
	Log& operator =(const Log&);
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
