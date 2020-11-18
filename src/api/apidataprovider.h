
#ifndef SRC_API_APIDATAPROVIDER_H_
#define SRC_API_APIDATAPROVIDER_H_

#include "feti4i.h"

#include <vector>
#include <functional>

namespace espreso {

class Kernel;

class APIDataProvider {
public:
	APIDataProvider();
	~APIDataProvider();

	int nodesSize();
	int matrixType();
	int DOFs();

	void prepare();
	void fillMatrix(std::function<void(FETI4IInt, FETI4IInt, FETI4IInt*, FETI4IReal*)> add);
	void fillDirichlet(std::vector<FETI4IInt> &indices, std::vector<FETI4IReal> &values);
	void fillL2G(std::vector<FETI4IInt> &l2g);
	void fillNeighbors(std::vector<FETI4IMPIInt> &neighbors);
	void fillRHS(std::vector<FETI4IReal> &rhs);

	void storeSolution(std::vector<FETI4IReal> &solution);

protected:
	Kernel *kernel;
	std::vector<double> rhs;
};

}

#endif /* SRC_API_APIDATAPROVIDER_H_ */
