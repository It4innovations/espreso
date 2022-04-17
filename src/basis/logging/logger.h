
#ifndef SRC_BASIS_LOGGING_LOGGER_H_
#define SRC_BASIS_LOGGING_LOGGER_H_

#include "verbosity.h"
#include <string>
#include <type_traits>

namespace espreso {

class LoggerBase {
public:
	int size = 0;
	VerboseArg** args;

	virtual void initOutput() =0;
	virtual void always() =0;
	virtual void start(const char* n, const char* section) =0;
	virtual void checkpoint(const char* n) =0;
	virtual void accumulated(const char* n) =0;
	virtual void end(const char* n) =0;
	virtual void ln() =0;
	virtual void nextLoadStep(int step) =0;
	virtual void param(const char* n, const int &value) =0;
	virtual void param(const char* n, const long &value) =0;
	virtual void param(const char* n, const long unsigned int &value) =0;
	virtual void param(const char* n, const double &value) =0;
	virtual void param(const char* n, const char* value) =0;
	virtual void finish() =0;
	virtual void output(const char* msg, VerboseArg::COLOR color) =0;
	virtual void error(const char* msg) =0;

	virtual ~LoggerBase() { delete[] args; }
protected:
	LoggerBase(int n): size(n), args(new VerboseArg*[n]) { }
};


#define __ES__PACK(...) __VA_ARGS__
#define __ES__FORALL(fnc, header, call)        \
	template<class Logger, class... Other>     \
	typename                                   \
	std::enable_if<sizeof...(Other)>::type     \
	fnc(header)                                \
	{                                          \
		fnc<Logger>(call);                     \
		fnc<Other...>(call);                   \
	}                                          \
	void fnc(header)                           \
	{                                          \
		fnc<Loggers...>(call);                 \
	}

template <class... Loggers>
class Logger: public LoggerBase, public Loggers... {

private:

	template<class Logger>
	void iterate()
	{
		args[size++] = &(static_cast<Logger&>(*this));
	}

	__ES__FORALL(iterate, , );

	template<class TSearch, class Logger>
	void search(VerboseArg* &logger)
	{
		if (std::is_same<TSearch, Logger>::value) {
			logger = static_cast<Logger*>(this);
		}
	}

	template<class TSearch, class Logger, class... Other>
	typename std::enable_if<sizeof...(Other)>::type
	search(VerboseArg* &logger)
	{
		search<TSearch, Logger>(logger);
		search<TSearch, Other...>(logger);
	}

	template<class TSearch>
	void search(VerboseArg* &logger)
	{
		search<TSearch, Loggers...>(logger);
	}

public:

	Logger(): LoggerBase(sizeof...(Loggers))
	{
		size = 0;
		iterate();
	}

	template <class TLogger>
	VerboseArg* get()
	{
		VerboseArg* instance = 0;
		search<TLogger, Loggers...>(instance);
		return instance;
	}

	template<class Logger>
	void initOutput()
	{
		Logger::initOutput();
	}

	template<class Logger>
	void always()
	{
		Logger::always = Logger::level;
	}

	template<class Logger>
	void start(const char* n, const char* section)
	{
		++Logger::level;
		if (Logger::isAllowed()) {
			Logger::start(n, section);
		}
	}

	template<class Logger>
	void checkpoint(const char* n)
	{
		if (Logger::isAllowed()) {
			Logger::checkpoint(n);
		}
	}

	template<class Logger>
	void accumulated(const char* n)
	{
		if (Logger::isAllowed()) {
			Logger::accumulated(n);
		}
	}

	template<class Logger>
	void end(const char* n)
	{
		if (Logger::isAllowed()) {
			Logger::end(n);
		}
		Logger::finishing = true;
	}

	template<class Logger>
	void param(const char* n, const int &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const long &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const long unsigned int &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const double &value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void param(const char* n, const char* value)
	{
		if (Logger::isAllowed()) {
			Logger::param(n, value);
		}
	}

	template<class Logger>
	void ln()
	{
		if (Logger::isAllowed()) {
			Logger::ln();
		}
		if (Logger::finishing) {
			Logger::finishing = false;
			--Logger::level;
		}
		if (Logger::always == Logger::level) {
			Logger::always = 0;
		}
	}

	template<class Logger>
	void nextLoadStep(int step)
	{
		Logger::nextLoadStep(step);
	}

	template<class Logger>
	void finish()
	{
		Logger::finish();
	}

	template<class Logger>
	void output(const char* msg, VerboseArg::COLOR color)
	{
		if (Logger::verbosity) {
			Logger::output(msg, color);
		}
	}

	template<class Logger>
	void error(const char* msg)
	{
		Logger::error(msg);
	}

	~Logger() {}

	__ES__FORALL(initOutput, , );
	__ES__FORALL(always, , );
	__ES__FORALL(start, __ES__PACK(const char* n, const char* section), __ES__PACK(n, section));
	__ES__FORALL(checkpoint, const char* n, n);
	__ES__FORALL(accumulated, const char* n, n);
	__ES__FORALL(end, const char* n, n);
	__ES__FORALL(ln, , );
	__ES__FORALL(nextLoadStep, int step, step);
	__ES__FORALL(param, __ES__PACK(const char* n, const int &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const long &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const long unsigned int &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const double &value), __ES__PACK(n, value));
	__ES__FORALL(param, __ES__PACK(const char* n, const char* value), __ES__PACK(n, value));
	__ES__FORALL(finish, , );
	__ES__FORALL(output, __ES__PACK(const char* msg, VerboseArg::COLOR color), __ES__PACK(msg, color));
	__ES__FORALL(error, const char* msg, msg);
};
}



#endif /* SRC_BASIS_LOGGING_LOGGER_H_ */
