
#ifndef BASIS_LOGGING_TIMEEVAL_H_
#define BASIS_LOGGING_TIMEEVAL_H_

#include <string>
#include <vector>

namespace espreso {

#ifdef FETI_DEBUG
#define DEBUGOUT std::cout
#else
#define DEBUGOUT if (true) ; else std::cout
#endif

struct TimeEvent
{
    friend class TimeEval;

    TimeEvent(std::string name);

    static double time();

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
    double getLastStat(double totalTime = 0.0);

    void printStatMPI(double totalTime = 0.0);
    void printLastStatMPI(double totalTime = 0.0);
    void printLastStatMPIPerNode(double totalTime = 0.0);

private:
    void evaluate();
    void evaluateMPI();

    std::string eventName;
    esint eventCount;
    std::vector<double> eventTime;

    esint name_length;
    esint val_length;

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
    void addPointerToEvent(TimeEvent *timeEvent);
    void printStats();
    void printStatsMPI();

    TimeEvent totalTime;
    TimeEvent remainingTime;

    std::string evalName;

    std::string eventName;
    std::vector<TimeEvent> timeEvents;
    std::vector<TimeEvent*> ptimeEvents;

};

}

#endif /* BASIS_LOGGING_TIMEEVAL_H_ */
