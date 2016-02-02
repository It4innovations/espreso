
#include "mkl.h"
#include "mpi.h"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <math.h>
#include <stack>
#include <ctime>

using std::vector;
using std::cout;
using std::string; 
using std::left;
using std::setw; 
using std::fixed;
using std::endl; 

#include "utils.h"

#pragma once

class TimeEvent
{
public:
	TimeEvent(string name); 
	TimeEvent(void); 
	~TimeEvent(void);

	string eventName; 
	eslocal eventCount;
	SEQ_VECTOR<double> eventTime; 

	eslocal name_length;
	eslocal val_length;

	double avgTime; 
	double sumTime; 
	double minTime; 
	double maxTime; 
	double stdDev; 

	double g_avgTime; 
	double g_minTime; 
	double g_maxTime; 
	double g_stdDev; 
	double g_sumTime;

	void SetName(string name);
	
	double getTime();
	void AddStart(double time); 
	void AddStart();
	void AddStartWithBarrier();
	void AddStartWOBarrier(double time);

	void AddEnd(double time); 
	void AddEnd(); 
	void AddEndWithBarrier();
	void AddEndWOBarrier(double time);

	void Reset(); 

	void Evaluate();
	void EvaluateMPI();
	
	void PrintStatMPI(); 
	void PrintStat(double totalTime); 
	void PrintStatMPI(double totalTime);
	void PrintLastStat(double totalTime);
	void PrintLastStatMPI(double totalTime);

	void PrintLastStatMPI_PerNode(double totalTime);


};


class TimeEval
{
public:

	TimeEval(void);
	TimeEval(string name);
	~TimeEval(void);

	TimeEvent totalTime; 
	TimeEvent remainingTime; 

	string evalName; 

	string eventName; 
	SEQ_VECTOR<TimeEvent> timeEvents; 

	void AddEvent(TimeEvent timeEvent); 
	void SetName(string name);
	void PrintStats();
	void PrintStatsMPI();
};

