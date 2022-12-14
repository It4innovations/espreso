
#ifndef SRC_BASIS_LOGGING_PROGRESSLOGGER_H_
#define SRC_BASIS_LOGGING_PROGRESSLOGGER_H_

#include "verbosity.h"
#include <cstdio>

namespace espreso {

#define BUFFER_SIZE 2048
#define PREBUFFER_SIZE 4 * BUFFER_SIZE

template <class TStream>
class ProgressLogger {

	char buffer[BUFFER_SIZE];

public:
	void start(const char* region, const char* section)
	{
		snprintf(buffer, BUFFER_SIZE, "%*s%s", static_cast<TStream*>(this)->level, " ", region);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void checkpoint(const char* region)
	{
		snprintf(buffer, BUFFER_SIZE, "%*s%s", static_cast<TStream*>(this)->level, " ", region);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void accumulated(const char* region)
	{
		checkpoint(region);
	}

	void end(const char* region)
	{
		snprintf(buffer, BUFFER_SIZE, "%*s%s", static_cast<TStream*>(this)->level, " ", region);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void param(const char* name, const int &value)
	{
		snprintf(buffer, BUFFER_SIZE, " [%s=%d]", name, value);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void param(const char* name, const long &value)
	{
		snprintf(buffer, BUFFER_SIZE, " [%s=%ld]", name, value);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void param(const char* name, const long unsigned int &value)
	{
		snprintf(buffer, BUFFER_SIZE, " [%s=%lu]", name, value);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void param(const char* name, const double &value)
	{
		snprintf(buffer, BUFFER_SIZE, " [%s=%f]", name, value);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void param(const char* name, const char* value)
	{
		snprintf(buffer, BUFFER_SIZE, " [%s=%s]", name, value);
		static_cast<TStream*>(this)->_output(buffer);
	}

	void ln()
	{
		static_cast<TStream*>(this)->_output("\n");
	}

	void nextLoadStep(int step)
	{
		// do nothing
	}

	void output(const char* msg, VerboseArg::COLOR color)
	{
		static_cast<TStream*>(this)->_color(msg, color);
	}

	void error(const char* msg)
	{
		static_cast<TStream*>(this)->_error(msg);
	}
};

class ProgressFileLogger: public ProgressLogger<ProgressFileLogger>, public Verbosity<ProgressFileLogger, 'v'> {
public:
	void initOutput()
	{
		setLogFile();
	}

	void finish()
	{

	}

	void _output(const char* msg)
	{
		if (rank == 0) {
			if (out) {
				fprintf(out, "%s", msg);
			} else {
				bsize += snprintf(prebuffer + bsize, PREBUFFER_SIZE - bsize, "%s", msg);
			}
		}
	}

	void _color(const char* msg, VerboseArg::COLOR color)
	{
		_output(msg);
	}

	void _error(const char* msg)
	{
		if (err) {
			fprintf(err, "%s", msg);
		}
	}

	ProgressFileLogger();
	~ProgressFileLogger();
	void setLogFile();
	void closeLogFile();

protected:
	FILE *out, *err;
	size_t bsize;
	char prebuffer[PREBUFFER_SIZE] = { 0 };
};

class ProgressTerminalLogger: public ProgressLogger<ProgressTerminalLogger >, public Verbosity<ProgressTerminalLogger, 'v'> {
public:
	void initOutput()
	{

	}

	void finish()
	{

	}

	void _output(const char* msg)
	{
		if (grank == 0) {
			printf("%s", msg);
			fflush(stdout);
		}
	}

	void _color(const char* msg, VerboseArg::COLOR color)
	{
		if (grank == 0) {
			if (color == VerboseArg::COLOR::WHITE) {
				_output(msg);
			} else {
				printf("%s%s\x1b[0m", getColor(color), msg);
				fflush(stdout);
			}
		}
	}

	void _error(const char* msg)
	{
		fprintf(stderr, "\x1b[31m%s\x1b[0m", msg);
		fflush(stderr);
	}

protected:
	const char* getColor(VerboseArg::COLOR color)
	{
		switch (color) {
		case VerboseArg::COLOR::RED:     return "\x1b[31m";
		case VerboseArg::COLOR::GREEN:   return "\x1b[32m";
		case VerboseArg::COLOR::YELLOW:  return "\x1b[33m";
		case VerboseArg::COLOR::BLUE:    return "\x1b[34m";
		case VerboseArg::COLOR::MAGENTA: return "\x1b[35m";
		case VerboseArg::COLOR::CYAN:    return "\x1b[36m";
		default:
			return "\x1b[0m";
		}
	}
};

}


#endif /* SRC_BASIS_LOGGING_PROGRESSLOGGER_H_ */
