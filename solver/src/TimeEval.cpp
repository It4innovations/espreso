#include "TimeEval.h"
#include "omp.h"

TimeEvent::TimeEvent(void)
{
	eventCount  = 0; 
	name_length = 50;
	val_length  = 15;

}


TimeEvent::TimeEvent(string name)
{
	eventName   = name;  
	eventCount  = 0; 
	name_length = 50;
	val_length  = 15;

}


TimeEvent::~TimeEvent(void)
{
}

void TimeEvent::SetName(string name)
{
	eventName = name; 
}

double TimeEvent::getTime() {
	return omp_get_wtime();
}

void TimeEvent::AddStart() {
#ifdef TM_BLOCK_START
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	double time = getTime();
	eventTime.push_back(time); 
} 

void TimeEvent::AddStart(double time) { 
#ifdef TM_BLOCK_START
	MPI_Barrier(MPI_COMM_WORLD);
	time = getTime();
#endif	
	eventTime.push_back(time); 
} 


void TimeEvent::AddStartWithBarrier() {
	MPI_Barrier(MPI_COMM_WORLD);
	double time = getTime();
	eventTime.push_back(time); 
} 

void TimeEvent::AddEndWithBarrier() {
	MPI_Barrier(MPI_COMM_WORLD);
	double time = getTime();
	eventTime[eventTime.size()-1] = time - eventTime[eventTime.size()-1]; 
	eventCount++; 
}


void TimeEvent::AddStartWOBarrier(double time) {
	eventTime.push_back(time); 
} 

void TimeEvent::AddEndWOBarrier(double time) {
	eventTime[eventTime.size()-1] = time - eventTime[eventTime.size()-1]; 
	eventCount++; 
}


void TimeEvent::AddEnd() {
#ifdef TM_BLOCK_END
	MPI_Barrier(MPI_COMM_WORLD);
#endif 
	double time = getTime();
	eventTime[eventTime.size()-1] = time - eventTime[eventTime.size()-1]; 
	eventCount++; 
}

void TimeEvent::AddEnd(double time) {
#ifdef TM_BLOCK_END
	MPI_Barrier(MPI_COMM_WORLD);
	time = getTime();
#endif
	eventTime[eventTime.size()-1] = time - eventTime[eventTime.size()-1]; 
	eventCount++; 
}

void TimeEvent::Reset() {
	eventCount = 0; 
	eventTime.clear();

	avgTime = 0.0;
	sumTime = 0.0; 
	minTime = 0.0; 
	maxTime = 0.0; 
	stdDev  = 0.0; 

	g_avgTime = 0.0; 
	g_minTime = 0.0; 
	g_maxTime = 0.0; 
	g_stdDev  = 0.0; 
	g_sumTime = 0.0;
}

void TimeEvent::Evaluate() {
	sumTime = 0; 
	avgTime = 0; 
	minTime = 1000000000000000000000000000.0; //pozor   
	maxTime = 0; 
	stdDev  = 0; 

	for (int i = 0; i < eventCount; i++) {
		double d_time = eventTime[i]; 
		sumTime += d_time; 
		if (d_time < minTime )
			minTime = d_time; 
		if (d_time > maxTime)
			maxTime = d_time;
	}
	
	avgTime = sumTime / eventCount; 

	double E=0;
	for (int i = 0; i < eventCount; i++)
		E+=(eventTime[i] - avgTime)*(eventTime[i] - avgTime);
	
	stdDev = sqrt(1/eventCount*E);
}

void TimeEvent::PrintStat(double totalTime) {
	Evaluate(); 

	cout << setw(name_length) << left << eventName << " avg.: " << fixed << setw(val_length) << avgTime << "min.: " << setw(val_length) << minTime << "max.: " << setw(val_length) << maxTime;

	if (totalTime == 0.0) {
		  ;
	} else {
		cout << "% of avg tot: " << setw(val_length) << 100.0 * avgTime/totalTime;
	}

	cout << endl; 

}


void TimeEvent::PrintLastStat(double totalTime) {

	avgTime = eventTime[eventTime.size()-1];
	
	cout << setw(name_length) << left << eventName << " avg.: " << fixed << setw(val_length) << avgTime << "min.: " << setw(val_length) << " -- " << "max.: " << setw(val_length) << " -- "; 

	if (totalTime == 0.0) {
		 ;
	} else {
		cout << "% of avg tot: " << setw(val_length) << 100.0 * avgTime/totalTime; 
	}
	
	cout << endl; 
}


void TimeEvent::EvaluateMPI() {
	Evaluate();

	int rank, size; 

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
	int mpi_root = 0;

	MPI_Reduce(&avgTime, &g_avgTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	g_avgTime = g_avgTime / size; 

	MPI_Reduce(&sumTime, &g_sumTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	g_sumTime = g_sumTime / size; 

	MPI_Reduce(&minTime, &g_minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); 
	MPI_Reduce(&maxTime, &g_maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 

}


void TimeEvent::PrintStatMPI(double totalTime) {
	EvaluateMPI();

	int rank; 
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */

	
	
	if(rank == 0) {
		cout << setw(name_length) << left << eventName << 
			" avg.: " << setw(val_length) << fixed << g_avgTime << 
			" min.: " << setw(val_length) << g_minTime << 
			" max.: " << setw(val_length) << g_maxTime << 
			" sum.: " << setw(val_length) << g_sumTime <<
			" count: " << setw(val_length) << eventCount;
		if (totalTime == 0.0) {
			  ;
		} else {
			cout << "% of avg tot: " << setw(val_length) << 100.0 * g_avgTime/totalTime;  
		}
		cout << endl;
	}
}


void TimeEvent::PrintLastStatMPI(double totalTime) {

	double d_time = eventTime[eventTime.size()-1]; 
	
	int rank, size; 

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
	int mpi_root = 0;

	MPI_Reduce(&d_time, &g_avgTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	g_avgTime= g_avgTime / size; 

	MPI_Reduce(&d_time, &g_minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); 
	MPI_Reduce(&d_time, &g_maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 

	if(rank == 0) {
		cout << setw(name_length) << left << eventName << 
			"last: " << fixed << setw(val_length) << g_avgTime << 
			"min.: " << setw(val_length) << g_minTime << 
			"max.: " << setw(val_length) << g_maxTime;
		if (totalTime == 0.0) {
			  ;
		} else {
			cout << "% of avg tot: " << setw(val_length) << 100.0 * g_avgTime/totalTime;  
		}
		cout << endl;
	}
}


void TimeEvent::PrintLastStatMPI_PerNode(double totalTime) {

	double d_time = eventTime[eventTime.size()-1]; 

	SEQ_VECTOR <double> d_all_times; 

	int rank, size; 

	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
	int mpi_root = 0;
	
	if(rank == 0)
		d_all_times.resize(size);
	else 
		d_all_times.resize(1);
		
	MPI_Gather(&d_time,1,MPI_DOUBLE,&d_all_times[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Reduce(&d_time, &g_avgTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	g_avgTime= g_avgTime / size; 

	MPI_Reduce(&d_time, &g_minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); 
	MPI_Reduce(&d_time, &g_maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); 

	if(rank == 0) {
		cout << setw(name_length) << left << eventName << 
			"last: " << fixed << setw(val_length) << g_avgTime << 
			"min.: " << setw(val_length) << g_minTime << 
			"max.: " << setw(val_length) << g_maxTime;
		if (totalTime == 0.0) {
			;
		} else {
			cout << "% of avg tot: " << setw(val_length) << 100.0 * g_avgTime/totalTime;  
		}
		cout << endl;

		for (int i = 0; i < size; i++) {
			cout << fixed << setw(3) << "R: " << setw(5) << i << setw(15) << d_all_times[i]; 
				
			if ((i+1) % 10 == 0)
				cout << endl; 
		}
		cout << endl;
	}
} 



////////////////////////////////////////////////////////////////////////////////////////////////

TimeEval::TimeEval(void) {
	totalTime     = TimeEvent(string("Total "));
	remainingTime = TimeEvent(string("Remaining "));
}


TimeEval::TimeEval(string name) {
	evalName = name; 
	totalTime =     TimeEvent(name + string("- Total "));
	remainingTime = TimeEvent(name + string("- Remaining "));

}

TimeEval::~TimeEval(void) {
	
}
 

void TimeEval::SetName(string name)
{
	evalName = name; 
}


void TimeEval::AddEvent(TimeEvent timeEvent) {
	timeEvents.push_back(timeEvent); 
}

void TimeEval::PrintStats() {
	totalTime.Evaluate(); 

	for (int i = 0; i < timeEvents.size(); i++) {
		timeEvents[i].PrintStat(totalTime.avgTime);
	}

	totalTime.PrintStat(totalTime.avgTime); 
}

void TimeEval::PrintStatsMPI() {
	int rank; 
	int size; 
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
	int mpi_root = 0;
	
	double sum_avg_time = 0;  
	
	if (rank == 0) {
		cout << endl;
		cout << "*****************************************************************************************************************************************************************************************" << endl; 
		cout << "*** " << evalName << " ***" << endl;
		cout << "*****************************************************************************************************************************************************************************************" << endl; 
	}
	totalTime.EvaluateMPI();

	for (int i = 0; i < timeEvents.size(); i++) {
		timeEvents[i].PrintStatMPI(totalTime.g_avgTime);
		sum_avg_time += timeEvents[i].avgTime; 
	}
	
	if (rank == 0) 
		if (rank == 0) 
			cout << "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

	totalTime.PrintStatMPI(totalTime.g_avgTime); 	
	//remainingTime.AddStart(sum_avg_time);
	//remainingTime.AddEnd  (totalTime.avgTime);
	//remainingTime.PrintStatMPI(totalTime.g_avgTime);
	if (rank == 0) 
		cout << "*****************************************************************************************************************************************************************************************" << endl << endl; 
}
