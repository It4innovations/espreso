/*! \file   elapsed_time.cpp
    \brief  time esurment functions
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 4th 2013
    \date   Jul. 12th 2015
    \date   Feb. 29th 2016
*/

#ifndef _elapsed_time_

#ifdef CLOCK_GETTIME

#  ifdef _MSC_VER   // added by Yann Collette
#  include <windows.h>

typedef struct timeval elapsed_t;

#define COPYTIME(a, b) ((a).tv_sec = (b).tv_sec);\
((a).tv_usec = (b).tv_usec)

#  else // ! _MSC_VER == Linux
#include <time.h>
typedef struct timespec elapsed_t;
#define COPYTIME(a, b) ((a).tv_sec = (b).tv_sec);\
((a).tv_nsec = (b).tv_nsec)

#endif

#else   /* #ifdef CLOCK_GETTIME */

#  ifdef GETRUSAGE
#    include <sys/time.h>
#    include <sys/resource.h>
typedef struct rusage elapsed_t;

#define COPYTIME(a, b) ((a).ru_utime.tv_sec = (b).ru_utime.tv_sec); \
((a).ru_utime.tv_usec = (b).ru_utime.tv_usec); \
((a).ru_stime.tv_sec = (b).ru_stime.tv_sec); \
((a).ru_stime.tv_usec = (b).ru_stime.tv_usec) 

#  else /* #ifdef GETTIMEOFDAY */

#ifdef CLOCK // for NEC SX-ACE
#include <time.h>
typedef clock_t elapsed_t;
#define COPYTIME(a, b) (a = b);

#else

#include <sys/time.h>
typedef struct timeval elapsed_t;

#define COPYTIME(a, b) ((a).tv_sec = (b).tv_sec); \
((a).tv_usec = (b).tv_usec)

# endif  /* #ifdef CLOCK */
#  endif /* #ifdef GETRUSAGE */
#endif /* #ifdef CLOCK_GETTIME */

void get_realtime(elapsed_t *tm);
double convert_time(elapsed_t time0, elapsed_t time1);
int convert_sec(elapsed_t t);
int convert_microsec(elapsed_t t);

#define _elapsed_time_
#endif

