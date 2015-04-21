
#ifndef SOLID45_H_
#define SOLID45_H_

#include "utility.h"

class CSolid45 {
public:
	  longInt inod_glob[8];
	  int indSubdomain;
	  int ordinalNumber;

public:
	  CSolid45();
	  virtual ~CSolid45();
};

#endif /* SOLID45_H_ */
