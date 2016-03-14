#include "../generic/utils.h"

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

eslocal SaveBinVectorDouble(SEQ_VECTOR <double> & Vector, string filename) {

	std::ofstream file (filename.c_str(), std::ios::out | std::ofstream::binary);

	if ( file.is_open() ) {

		file.write((char*) &Vector[0], Vector.size() * sizeof(double) / sizeof(char)); 
		file.close();

		return 0; 

	} else {

		cout << "File " << filename << " not found ! " << endl; 

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



eslocal LoadBinVectorInt(SEQ_VECTOR <eslocal> & Vector, string filename) {

	// open the file:
	std::streampos fileSize;
	std::ifstream file (filename.c_str(), std::ios::binary);

	if ( file.is_open() ) {
		// get its size:
		file.seekg(0, std::ios::end);
		fileSize = file.tellg();
		file.seekg(0, std::ios::beg);

		// read the data:
		Vector.resize(fileSize / sizeof(eslocal));

		file.read((char*) &Vector[0], fileSize);

		file.close();

		return 0; 
	} else {
		cout << "File " << filename << " not found ! " << endl; 
		return -1; 
	}
}

eslocal LoadBinVecVec(SEQ_VECTOR <SEQ_VECTOR <eslocal> > & outputVecVec, string filename) {

	std::ifstream in (filename.c_str(), std::ios::binary);
	char delim = ';'; 
	string line, field;

	if ( in.is_open() ) {

		// Get parameters 
		getline(in,line);
		std::stringstream paramss(line);

		getline(paramss,field,delim); 
		eslocal rows = atoi(field.c_str());		// get num of rows

		getline(paramss,field,delim); 
		eslocal cols = atoi(field.c_str());		// get num of columns

		// Get data 

		for (eslocal r_i = 0; r_i < rows; r_i++) {
			SEQ_VECTOR <eslocal> tmp_coor (cols);
			in.read((char*) &tmp_coor [0], cols*sizeof(eslocal));
			outputVecVec.push_back(tmp_coor); 
		}

		in.close(); 

		return 0; 
	} else {
		cout << "File " << filename << " not found ! " << endl; 
		return -1; 
	}
}

eslocal LoadBinVecVec(SEQ_VECTOR <SEQ_VECTOR <double> > & outputVecVec, string filename) {

	std::ifstream in (filename.c_str(), std::ios::binary);
	char delim = ';'; 
	string line, field;

	if ( in.is_open() ) {

		// Get parameters 
		getline(in,line);
		std::stringstream paramss(line);

		getline(paramss,field,delim); 
		eslocal rows = atoi(field.c_str());		// get num of rows

		getline(paramss,field,delim); 
		eslocal cols = atoi(field.c_str());		// get num of columns

		// Get data 

		for (eslocal r_i = 0; r_i < rows; r_i++) {
			SEQ_VECTOR <double> tmp_coor (cols); 
			in.read((char*) &tmp_coor [0], cols*sizeof(double));
			outputVecVec.push_back(tmp_coor); 
		}

		in.close(); 

		return 0; 
	} else {
		cout << "File " << filename << " not found ! " << endl; 
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
		cout << endl << "Thread " << omp_get_thread_num() << " - Printing vector : " << name << endl; 
		for (eslocal i = 0; i < vec.size(); i++) {
			cout << vec[i] << endl; 
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
		cout << endl << "Thread " << omp_get_thread_num() << " - Printing vector : " << name << endl; 
		for (eslocal i = 0; i < vec.size(); i++) {
			cout << vec[i] << endl; 
		}
	}
}

static eslocal parseLine_u(char* line){
	eslocal i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}

void GetProcessMemoryStat_u ( ) {

#ifndef WIN32

	FILE* file = fopen("/proc/self/status", "r");
	eslocal result = -1;
	char line[128];


	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmRSS:", 6) == 0){
			result = parseLine_u(line);
			break;
		}
	}
	fclose(file);

	int MPIrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

	//cout << endl; 
	//cout << " ******************************************************************************************************************************* " << endl; 
	cout << " - Memory used by process " << MPIrank << " : " << result / 1024.0 << " MB"<< endl; 
	//cout << " ******************************************************************************************************************************* " << endl; 
	//cout << endl; 

#endif

}

void GetMemoryStat_u( ) 
{
#ifndef WIN32

	struct sysinfo memInfo;
	sysinfo (&memInfo);

	long long totalPhysMem;
	long long physMemUsed; 

	totalPhysMem = memInfo.totalram;
	//Multiply in next statement to avoid int overflow on right hand side...
	totalPhysMem *= memInfo.mem_unit;

	physMemUsed	= memInfo.totalram - memInfo.freeram;
	//Multiply in next statement to avoid int overflow on right hand side...
	physMemUsed *= memInfo.mem_unit;

	//	cout << endl; 
	//	cout << " ******************************************************************************************************************************* " << endl; 
	//	cout << " *** Memory Info ... " << endl; 
	//	cout << "  - Total RAM memory : " << totalPhysMem << endl; 
	//	cout << "  - Used RAM  memory : " << physMemUsed << endl;
	//	cout << "  - Usage            : " << 100.0 * (double)physMemUsed/(double)totalPhysMem<< " % " << endl ; 
	//	cout << " ******************************************************************************************************************************* " << endl; 
	//	cout << endl; 
	cout << " - Total used RAM : " << 100.0 * (double)physMemUsed/(double)totalPhysMem<< " %  - " << physMemUsed/1024/1024 << " MB of " << totalPhysMem/1024/1024 << " MB" << endl;

#endif
}

double GetProcessMemory_u ( ) {

#ifndef WIN32

	FILE* file = fopen("/proc/self/status", "r");
	eslocal result = -1;
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


}
// **** END - Uncategorized functions ********************************
// *******************************************************************
