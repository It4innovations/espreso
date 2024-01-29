
#ifndef SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_
#define SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_

#include "feti/feti.h"

namespace espreso {

template <typename T>
struct Regularization {

	Regularization(FETI<T> &feti): feti(feti) {}
	virtual ~Regularization() {}

	void set(step::Step &step);
	void update(step::Step &step);

protected:
	void setAlgebraic();
	void updateAlgebraic();

	virtual void setAnalytic() { setAlgebraic(); }
	virtual void updateAnalytic() { updateAlgebraic(); }

	void getFixPoints(std::vector<esint> &fixPoints, int domain);

	FETI<T> &feti;
};

}

#endif /* SRC_ANALYSIS_LINEARSYSTEM_REGULARIZATION_REGULARIZATION_H_ */
