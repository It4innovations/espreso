
#ifndef SRC_OUTPUT2_RESULTSTORE_H_
#define SRC_OUTPUT2_RESULTSTORE_H_

#include <string>
#include <map>

#include "store.h"
#include "../basis/point/point.h"

namespace espreso {

class Region;
struct Solution;

namespace output {

struct DataArrays {
	std::map<std::string, std::vector<eslocal>* > pointDataInteger, elementDataInteger;
	std::map<std::string, std::vector<double>* > pointDataDouble, elementDataDouble;

	~DataArrays();
};

class ResultStore: public Store {

public:
	virtual void storeSettings(const Step &step);
	virtual void storeSettings(size_t steps);
	virtual void storeSettings(const std::vector<size_t> &steps);

	virtual void storeSolution(const Step &step, const Solution &solution);

	virtual ~ResultStore() {};

protected:
	ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);

	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data) =0;
	virtual void compose(const std::string &name, const std::vector<std::string> &names) {};

	virtual void preprocessing();
	virtual void regionPreprocessing(const espreso::Region &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements);
	virtual void regionData(size_t step, const espreso::Region &region, DataArrays &data);

	virtual void coordinatePreprocessing(const std::vector<std::vector<eslocal> > &indices, std::vector<double> &coordinates, std::vector<size_t> &offsets);

	std::vector<double> _coordinates; // x1, y1, z1, x2, y2, z2, ...
	std::vector<eslocal> _elementsTypes;  // code1, code2, ...
	std::vector<eslocal> _elementsNodes;  // nodes1, nodes2, ...
	std::vector<eslocal> _elements;  // n11, n12, n13, ..., n21, n22, n23, ...

	Point _clusterCenter;
	std::vector<Point> _domainsCenters;
};

}
}

#endif /* SRC_OUTPUT2_RESULTSTORE_H_ */
