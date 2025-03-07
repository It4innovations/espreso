
#ifndef SRC_ESINFO_ESLOG_HPP_
#define SRC_ESINFO_ESLOG_HPP_

#include "eslog.h"
#include <cstdio>

namespace espreso {
namespace eslog {

#define BUFFER_SIZE 2048
extern char buffer[BUFFER_SIZE];

template <typename... Args>
void info(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		info(buffer);
	}
}

template <typename... Args>
void solver(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		solver(buffer);
	}
}

template <typename... Args>
void linearsolver(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		linearsolver(buffer);
	}
}

template <typename... Args>
void duration(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		duration(buffer);
	}
}

template <typename... Args>
void warning(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		warning(buffer);
	}
}

template <typename... Args>
void storedata(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		storedata(buffer);
	}
}

template <typename... Args>
[[noreturn]] void failure(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		failure(buffer);
	}
}

template <typename... Args>
[[noreturn]] void internalFailure(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		internalFailure(buffer);
	}
}

template <typename... Args>
[[noreturn]] void error(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		error(buffer);
	}
}

template <typename... Args>
[[noreturn]] void globalerror(const char* format, Args... args)
{
	#pragma omp critical(eslogbuffer)
	{
		if (BUFFER_SIZE < snprintf(buffer, BUFFER_SIZE, format, args...)) {
			buffer[BUFFER_SIZE - 2] = '\n';
		}
		globalerror(buffer);
	}
}

}
}



#endif /* SRC_ESINFO_ESLOG_HPP_ */
