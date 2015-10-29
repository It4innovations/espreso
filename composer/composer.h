
#ifndef COMPOSER_COMPOSER_H_
#define COMPOSER_COMPOSER_H_

#include "essolver.h"

namespace composer {

class Composer {

public:
	virtual void element(SparseMatrix &K, size_t subdomain) = 0;
	virtual void subdomainBoundary() = 0;
	virtual void clusterBoundary() = 0;

	virtual ~Composer() {};
protected:
	Composer() {};
};

}


#endif /* COMPOSER_COMPOSER_H_ */
