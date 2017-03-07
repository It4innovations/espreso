
#ifndef SRC_OUTPUT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULTSTORE_H_

#include <string>
#include <map>

#include "../basis/point/point.h"
#include "../assembler/step.h"

#include "store.h"

namespace espreso {

class Element;
class Mesh;
class Region;
enum class ElementType;

namespace output {

struct DataArrays {
	std::map<std::string, std::pair<size_t, std::vector<eslocal>* > > pointDataInteger, elementDataInteger;
	std::map<std::string, std::pair<size_t, std::vector<double>* > > pointDataDouble, elementDataDouble;

	~DataArrays();
};

class ResultStore: public Store {

public:
	const OutputConfiguration& configuration() const { return _configuration; }

	// TODO: remove it after removing old physics
	void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType);

	virtual void storeSettings(const Step &step);
	virtual void storeSettings(size_t steps);
	virtual void storeSettings(const std::vector<size_t> &steps);

	virtual void storeFETIData(const Step &step, const Instance &instance);

	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution);
	virtual void finalize();

	virtual ~ResultStore();

protected:
	ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path);
	virtual void preprocessing();

	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const DataArrays &data) =0;
	virtual void store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution) =0;

	virtual void linkClusters(const std::string &root, const std::string &name, const DataArrays &data) =0;
	virtual void linkClusters(const std::string &root, const std::string &name, const std::vector<Solution*> &solution) =0;

	virtual void linkSteps(const std::string &name, const std::vector<std::pair<std::string, Step> > &steps) =0;

	const Mesh *_mesh;
	std::string _path;

	std::vector<double> _coordinates; // x1, y1, z1, x2, y2, z2, ...
	std::vector<eslocal> _elementsTypes;  // code1, code2, ...
	std::vector<eslocal> _elementsNodes;  // nodes1, nodes2, ...
	std::vector<eslocal> _elements;  // n11, n12, n13, ..., n21, n22, n23, ...

	Point _clusterCenter;
	std::vector<Point> _domainsCenters;

	std::vector<std::pair<std::string, Step> > _steps;

private:
	template <class TData>
	std::string store(const std::string &name, const Step &step, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const TData &data);

	void elementsPreprocessing(const std::vector<Element*> &region, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements);
	void regionData(size_t step, const espreso::Region &region, DataArrays &data);

	void coordinatePreprocessing(const std::vector<std::vector<eslocal> > &indices, std::vector<double> &coordinates, std::vector<size_t> &offsets);

	void storeFixPoints(const Step &step);
	void storeCorners(const Step &step);
	void storeDirichlet(const Step &step, const Instance &instance);
	void storeLambdas(const Step &step, const Instance &instance);


};

}
}

#endif /* SRC_OUTPUT_RESULTSTORE_H_ */
