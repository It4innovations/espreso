/*! \file   OptionLibrary.hpp
    \brief  compatibility for Microsoft compiler
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Feb. 23th 2013
    \date   Feb. 29th 2016
*/

#ifndef _COMPILER_OPTIONLIBRARY_H
# define _COMPILER_OPTIONLIBRARY_H

#ifdef _MSC_VER
#  include <process.h>
#else
#  include <stdlib.h>
#  include <unistd.h>
#endif

#ifdef SX_ACE
#include <string>
#include <cstdlib>
#endif
static inline double random_normalized()
{
#ifdef _MSC_VER
  return ((double)rand() / (double)RAND_MAX);
#else
  #ifdef SX_ACE
  return ((double)rand() / (double)RAND_MAX);
  #else
  return ((double)random() / (double)RAND_MAX);
  #endif
#endif
}

static inline bool random_bool()
{
#ifdef _MSC_VER
  double r = (double)rand() / (double)RAND_MAX;
#else
#ifdef SX_ACE
  double r = (double)rand() / (double)RAND_MAX;
#else
  double r = (double)random() / (double)RAND_MAX;
#endif
#endif
  return (r < 0.5 ? true : false);
}

static inline int get_process_id()
{
#ifdef _MSC_VER
  return (int)_getpid();
#else
  return (int)getpid();
#endif
}

// Intel compiler + older GNU C++ library may not have to_string()
// SX_ACE does not have
#ifdef NO_TO_STRING
inline std::string to_string(int num)
{
  char buf[256];
  sprintf(buf, "%d", num);
  std::string st = buf;
  return st;
}
#endif
#endif
