
#ifndef SRC_ESINFO_ESLOG_H_
#define SRC_ESINFO_ESLOG_H_

namespace espreso {

struct VerboseArg;
class LoggerBase;

namespace eslog {

void init(LoggerBase *logger);
void initFiles();
void reinit();
void printRunInfo(int *argc, char ***argv);
void finish();

void always();
void start(const char* name, const char* section);
void checkpoint(const char* name);
void accumulated(const char* name);
void end(const char* name);
void ln();
void nextStep(int step);

void startln(const char* name, const char* section);
void checkpointln(const char* name);
void accumulatedln(const char* name);
void endln(const char* name);

void param(const char* name, const int &value);
void param(const char* name, const long &value);
void param(const char* name, const long unsigned int &value);
void param(const char* name, const double &value);
void param(const char* name, const char* value);

void info(const char* msg);
void solver(const char* msg);
void linearsolver(const char* msg);
void duration(const char* msg);
void warning(const char* msg);
void storedata(const char* msg);
[[noreturn]] void failure(const char* msg);
[[noreturn]] void internalFailure(const char* msg);
[[noreturn]] void error(const char* msg);
[[noreturn]] void globalerror(const char* msg);

double time();
double duration();

extern LoggerBase *logger;

}
}



#endif /* SRC_ESINFO_ESLOG_H_ */
