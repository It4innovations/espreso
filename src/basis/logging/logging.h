
#ifndef BASIS_LOGGING_LOGGING_H_
#define BASIS_LOGGING_LOGGING_H_

#include <string>
#include <sstream>
#include <vector>

#define ESTEST(EVENT) if (!espreso::Test::report(EVENT))    ; else espreso::Test(EVENT).get()
#define ESINFO(EVENT) if (!espreso::Info::report(EVENT))    ; else espreso::Info(EVENT).get()
#define ESLOG(EVENT)  if (!espreso::Measure::report(EVENT)) ; else espreso::Measure(EVENT).get()

namespace espreso {

enum ESPRESOTest {
	TEST_FAILED,
	TEST_PASSED
};

enum TestEvent {
	MANDATORY,
	EVALUATION,
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
	GLOBAL_ERROR,
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
	MEASURE_LEVEL3,

	CLUSTER
};

struct OutputConfiguration;
struct Step;

class Test
{
public:
	Test& operator<<(const ESPRESOTest &test)
	{
		if (test == TEST_FAILED) { error = true; }
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
		switch (outputLevel) {
		case 0: return event < TEST_LEVEL0;
		case 1: return event < TEST_LEVEL1;
		case 2: return event < TEST_LEVEL2;
		case 3: return event < TEST_LEVEL3;
		default : return true;
		}
	};

	static void setLevel(size_t level) { outputLevel = level; }

protected:
	static size_t outputLevel;

	std::ostringstream os;
	bool error;
};

class Info
{
public:
	enum class TextColor {
		WHITE,
		RED,
		GREEN,
		YELLOW,
		BLUE
	};

	Info(InfoEvent event): event(event), color(TextColor::WHITE), _plain(false)
	{
		if (testLevel) {
			switch (event) {
			case CONVERGENCE:
				os << "CONVERGENCE: ";
				break;
			default:
				break;
			}
		}
	};
	~Info();

	enum InfoMode { FORMATTED, PLAIN };
	Info& operator<<(const InfoMode& mode)
	{
		_plain = mode == PLAIN;
		return *this;
	}
	Info& operator<<(const TextColor& color)
	{
		this->color = color;
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
	template<typename Tvalue>
	static std::string averageValues(const std::vector<Tvalue> &values);

	static bool report(InfoEvent event) {
		switch (outputLevel) {
		case 0: return event < VERBOSE_LEVEL0;
		case 1: return event < VERBOSE_LEVEL1;
		case 2: return event < VERBOSE_LEVEL2;
		case 3: return event < VERBOSE_LEVEL3;
		default : return true;
		}
	};

	static void setLevel(size_t level, size_t testing) { outputLevel = level; testLevel = testing; }

protected:
	static size_t outputLevel;
	static size_t testLevel;

	std::ostringstream os;
	InfoEvent event;
	TextColor color;
	bool _plain;
};

struct Checkpoint;

class Measure
{
public:
	static double time();

	Measure(MeasureEvent event): event(event) {};
	~Measure();

	static double processMemory();
	static double usedRAM();
	static double availableRAM();

	std::ostringstream& get() { return os; };

	static bool report(MeasureEvent event) {
		switch (outputLevel) {
		case 0: return event < MEASURE_LEVEL0;
		case 1: return event < MEASURE_LEVEL1;
		case 2: return event < MEASURE_LEVEL2;
		case 3: return event < MEASURE_LEVEL3;
		default : return true;
		}
	};

	static void setLevel(size_t level) { outputLevel = level; }

protected:
	static size_t outputLevel;

	static std::vector<Checkpoint> checkpoints;

	void evaluateCheckpoints();

	std::ostringstream os;
	MeasureEvent event;
};

class Logging {

public:
	static std::string outputRoot();
	static std::string prepareFile(const std::string &name);
	static std::string prepareFile(size_t subdomain, const std::string &name);

	static std::string path;
	static std::string name;
	static time_t time;
	static std::string debug;
	static std::ofstream log;
	static Step *step;
	static int rank;
};
}


#endif /* BASIS_LOGGING_LOGGING_H_ */
