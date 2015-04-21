
#ifndef UTILITY_
#define UTILITY_

#include <mpi.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <time.h>
#include <list>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>  // maybe better place ???
#include <vector>
#ifndef WIN32
#include <unistd.h>
#include "sys/types.h"
#include "sys/sysinfo.h"
#endif

#ifdef FLLOP_ENABLED
#include <fllopaif.h>
#define LogEventRegister(ename,eid) FllopAIFLogEventRegister(ename,eid)
#define LogEventBegin(eid)          FllopAIFLogEventBegin(eid)
#define LogEventEnd(eid)            FllopAIFLogEventEnd(eid)
extern int assemble_stage, solve_stage;
#else
#define LogEventRegister(ename,eid)
#define LogEventBegin(eid)
#define LogEventEnd(eid)
#endif

#define CONST_PI 3.14159265359

//
using namespace std;
#include "fortran_wrapper.h"
#ifdef flag_Fortran
#include "fortran.h"
#endif


//typedef int64_t longInt;
//typedef long long longInt;

typedef int longInt;
typedef int shortInt; 




#endif /* UTILITY_ */
