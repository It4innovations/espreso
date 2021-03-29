#include "utils.h"

#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"
#include "basis/utilities/utils.h"

#include <sstream>
#include <fstream>
#include <cstring>

#include "mkl.h"

#include <limits>
#include <cmath>

using std::endl; 


namespace espreso {
// *******************************************************************
// **** Uncategorized functions **************************************

//void LoadVectorInt(vector <int> & Vector, string filename) {
//
//	ifstream in (filename.c_str());
//	char delim = ';'; 
//	string line, field;
//
//	// Get line with values - INT 
//	getline(in,line);
//
//	stringstream ssI(line);
//	while (getline(ssI,field,delim))  // break line into comma delimited fields
//	{
//		Vector.push_back(atoi(field.c_str()));  // add each field to the 1D array
//	}
//
//	in.close(); 
//}

//void LoadVectorDouble(vector <double> & Vector, string filename) {
//
//	ifstream in (filename.c_str());
//	char delim = ';'; 
//	string line, field;
//
//	// Get line with values - INT 
//	getline(in,line);
//
//	stringstream ssI(line);
//	while (getline(ssI,field,delim))  // break line into comma delimitted fields
//	{
//		Vector.push_back(atof(field.c_str()));  // add each field to the 1D array
//	}
//
//	in.close();
//}

esint SaveBinVectorDouble(SEQ_VECTOR <double> & Vector, string filename) {

	std::ofstream file (filename.c_str(), std::ios::out | std::ofstream::binary);

	if ( file.is_open() ) {

		file.write((char*) &Vector[0], Vector.size() * sizeof(double) / sizeof(char)); 
		file.close();

		return 0; 

	} else {
		eslog::error("File not found.\n");


		return -1; 

	}

}

//void LoadBinVectorDouble(vector <double> & Vector, string filename) {
//	
//	// open the file:
//	streampos fileSize;
//	ifstream file (filename.c_str(), std::ios::binary);
//
//	// get its size:
//	file.seekg(0, std::ios::end);
//	fileSize = file.tellg();
//	file.seekg(0, std::ios::beg);
//
//	// read the data:
//	Vector.resize(fileSize / sizeof(double)); 
//	
//	file.read((char*) &Vector[0], fileSize);
//
//	file.close();
//}



esint LoadBinVectorInt(SEQ_VECTOR <esint> & Vector, string filename) {

	// open the file:
	std::streampos fileSize;
	std::ifstream file (filename.c_str(), std::ios::binary);

	if ( file.is_open() ) {
		// get its size:
		file.seekg(0, std::ios::end);
		fileSize = file.tellg();
		file.seekg(0, std::ios::beg);

		// read the data:
		Vector.resize(fileSize / sizeof(esint));

		file.read((char*) &Vector[0], fileSize);

		file.close();

		return 0; 
	} else {
		eslog::error("File not found.\n");
		return -1; 
	}
}

esint LoadBinVecVec(SEQ_VECTOR <SEQ_VECTOR <esint> > & outputVecVec, string filename) {

	std::ifstream in (filename.c_str(), std::ios::binary);
	char delim = ';'; 
	string line, field;

	if ( in.is_open() ) {

		// Get parameters 
		getline(in,line);
		std::stringstream paramss(line);

		getline(paramss,field,delim); 
		esint rows = atoi(field.c_str());		// get num of rows

		getline(paramss,field,delim); 
		esint cols = atoi(field.c_str());		// get num of columns

		// Get data 

		for (esint r_i = 0; r_i < rows; r_i++) {
			SEQ_VECTOR <esint> tmp_coor (cols);
			in.read((char*) &tmp_coor [0], cols*sizeof(esint));
			outputVecVec.push_back(tmp_coor); 
		}

		in.close(); 

		return 0; 
	} else {
		eslog::error("File not found.\n");
		return -1; 
	}
}

esint LoadBinVecVec(SEQ_VECTOR <SEQ_VECTOR <double> > & outputVecVec, string filename) {

	std::ifstream in (filename.c_str(), std::ios::binary);
	char delim = ';'; 
	string line, field;

	if ( in.is_open() ) {

		// Get parameters 
		getline(in,line);
		std::stringstream paramss(line);

		getline(paramss,field,delim); 
		esint rows = atoi(field.c_str());		// get num of rows

		getline(paramss,field,delim); 
		esint cols = atoi(field.c_str());		// get num of columns

		// Get data 

		for (esint r_i = 0; r_i < rows; r_i++) {
			SEQ_VECTOR <double> tmp_coor (cols); 
			in.read((char*) &tmp_coor [0], cols*sizeof(double));
			outputVecVec.push_back(tmp_coor); 
		}

		in.close(); 

		return 0; 
	} else {
		eslog::error("File not found.\n");
		return -1; 
	}
}






template <typename T>
void PrintVec(SEQ_VECTOR <T> vec, string name) {
#if DEBUG == 1
#ifdef _OPENMP
#pragma omp critical 
#endif
	{
		//ESINFO(ALWAYS_ON_ROOT) << "Thread " << omp_get_thread_num() << " - Printing vector : " << name;
		for (esint i = 0; i < vec.size(); i++) {
			//ESINFO(ALWAYS_ON_ROOT) << vec[i];
		}
	}
#endif // DEBUG
}

template <typename T>
void PrintVecND(SEQ_VECTOR <T> vec, string name) {
#ifdef _OPENMP
#pragma omp critical 
#endif
	{
//		//ESINFO(ALWAYS_ON_ROOT) << "Thread " << omp_get_thread_num() << " - Printing vector : " << name;
//		for (esint i = 0; i < vec.size(); i++) {
//			//ESINFO(ALWAYS_ON_ROOT) << vec[i];
//		}
	}
}

static esint parseLine_u(char* line){
	esint i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}

void GetProcessMemoryStat_u ( ) {

#ifndef WIN32
//
//	FILE* file = fopen("/proc/self/status", "r");
//	esint result = -1;
//	char line[128];
//
//
//	while (fgets(line, 128, file) != NULL){
//		if (strncmp(line, "VmRSS:", 6) == 0){
//			result = parseLine_u(line);
//			break;
//		}
//	}
//	fclose(file);

//	ESLOG(MEMORY) << " - Memory used by process " << info::mpi::rank << " : " << result / 1024.0 << " MB";


#endif

}

void GetMemoryStat_u( ) 
{
#ifndef WIN32

//	struct sysinfo memInfo;
//	sysinfo (&memInfo);
//
//	long long totalPhysMem;
//	long long physMemUsed;
//
//	totalPhysMem = memInfo.totalram;
//	//Multiply in next statement to avoid int overflow on right hand side...
//	totalPhysMem *= memInfo.mem_unit;
//
//	physMemUsed	= memInfo.totalram - memInfo.freeram;
//	//Multiply in next statement to avoid int overflow on right hand side...
//	physMemUsed *= memInfo.mem_unit;

//	ESLOG(MEMORY) << " - Total used RAM : " << 100.0 * (double)physMemUsed/(double)totalPhysMem<< " %  - " << physMemUsed/1024/1024 << " MB of " << totalPhysMem/1024/1024 << " MB";
#endif
}

double GetProcessMemory_u ( ) {

#ifndef WIN32

	FILE* file = fopen("/proc/self/status", "r");
	esint result = -1;
	char line[128];


	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmRSS:", 6) == 0){
			result = parseLine_u(line);
			break;
		}
	}
	fclose(file);

	return result / 1024.0;

#else 
	return 0.0;
#endif

}

int CompareVectors(SEQ_VECTOR <double> & vec_a, SEQ_VECTOR <double> & vec_b) {
	int err = 0;
//	double min_epsilon = std::numeric_limits<double>::epsilon();
	double epsilon = 1e-18;
	double diff = 0.0;

	for (size_t i = 0; i < vec_a.size(); i++)
	{
		diff = std::abs(vec_a[i] - vec_b[i]);
		// std::cout << diff << std::endl;
		if(diff >= epsilon) {
			err++;
		}
	}
	
	return err;
}

}
// **** END - Uncategorized functions ********************************
// *******************************************************************
