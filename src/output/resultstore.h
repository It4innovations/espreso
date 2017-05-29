
#ifndef SRC_OUTPUT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULTSTORE_H_

#include <string>
#include <map>

#include "../basis/point/point.h"
#include "../assembler/step.h"
#include "store.h"
#include "meshinfo.h"

namespace espreso {

class Element;
class Mesh;
class Region;
enum class ElementType;

namespace output {

struct RegionData;

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

	virtual std::string store(const std::string &name, const RegionData &regionData) =0;
	virtual std::string linkClusters(const std::string &root, const std::string &name, const RegionData &regionData) =0;
	virtual void linkSteps(const std::string &name, const std::vector<std::pair<Step, std::vector<std::string> > > &steps) =0;

protected:
	ResultStore(const OutputConfiguration &output, const Mesh *mesh, const std::string &path, MeshInfo::InfoMode mode = MeshInfo::EMPTY);

	const Mesh *_mesh;
	std::string _path;

	MeshInfo::InfoMode _mode;
	MeshInfo* _meshInfo;
	std::vector<std::pair<Step, std::vector<std::string> > > _solutions;
	std::vector<std::pair<Step, std::vector<std::string> > > _settings;
	std::vector<std::pair<Step, std::vector<std::string> > > _FETIdata;

private:
	std::vector<std::string> store(const std::string &name, const Step &step, const MeshInfo *meshInfo);

	void storeElementInfo(const Step &step);
	void storeFixPoints(const Step &step);
	void storeCorners(const Step &step);
	void storeDirichlet(const Step &step, const Instance &instance);
	void storeLambdas(const Step &step, const Instance &instance);


};

}
}

#endif /* SRC_OUTPUT_RESULTSTORE_H_ */
