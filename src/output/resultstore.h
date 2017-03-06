
#ifndef SRC_OUTPUT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULTSTORE_H_

#include <string>
#include <map>
#include <vector>

#include "../basis/point/point.h"

namespace espreso {

class Mesh;
class Region;
struct Step;
struct Solution;
struct OutputConfiguration;
enum class ElementType;

namespace output {

struct DataArrays {
	std::map<std::string, std::vector<eslocal>* > pointDataInteger, elementDataInteger;
	std::map<std::string, std::vector<double>* > pointDataDouble, elementDataDouble;

	~DataArrays();
};

class Store {

public:
	const OutputConfiguration& configuration() const { return _configuration; }

	virtual void storeSettings(const Step &step) =0;
	virtual void storeSettings(size_t steps) =0;
	virtual void storeSettings(const std::vector<size_t> &steps) =0;

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution) =0;
	virtual void finalize() =0;

	virtual ~Store() {};

protected:
	Store(const OutputConfiguration &configuration): _configuration(configuration) {};

	const OutputConfiguration &_configuration;
};

class ResultStore: public Store {

public:
	const OutputConfiguration& configuration() const { return _configuration; }

	// TODO: remove it after removing old physics
	void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType);

	virtual void storeSettings(const Step &step);
	virtual void storeSettings(size_t steps);
	virtual void storeSettings(const std::vector<size_t> &steps);

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution);

	virtual ~ResultStore() {};

protected:
	ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual void preprocessing();

	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data) =0;
	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution) =0;

	virtual void linkClusters(const std::string &root, const std::string &name, const DataArrays &data) =0;
	virtual void linkClusters(const std::string &root, const std::string &name, const std::vector<Solution*> &solution, size_t points, size_t cells) =0;

	virtual void linkSteps(const std::string &root, const std::string &name, const DataArrays &data) =0;
	virtual void linkSteps(const std::string &root, const std::string &name, const std::vector<Solution*> &solution) =0;

	const Mesh *_mesh;
	std::string _path;

	std::vector<double> _coordinates; // x1, y1, z1, x2, y2, z2, ...
	std::vector<eslocal> _elementsTypes;  // code1, code2, ...
	std::vector<eslocal> _elementsNodes;  // nodes1, nodes2, ...
	std::vector<eslocal> _elements;  // n11, n12, n13, ..., n21, n22, n23, ...

	Point _clusterCenter;
	std::vector<Point> _domainsCenters;

private:
	virtual void regionPreprocessing(const espreso::Region &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements);
	virtual void regionData(size_t step, const espreso::Region &region, DataArrays &data);

	virtual void coordinatePreprocessing(const std::vector<std::vector<eslocal> > &indices, std::vector<double> &coordinates, std::vector<size_t> &offsets);
};

}
}

#endif /* SRC_OUTPUT_RESULTSTORE_H_ */
