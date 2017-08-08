
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

struct RegionData;

class ResultStore: public Store {

	friend class AsyncStore;
	friend class AsyncStoreExecutor;

public:
	const OutputConfiguration& configuration() const { return _configuration; }

	virtual void storeSettings(size_t steps);
	virtual void storeFETIData(const Step &step, const Instance &instance);
	virtual void storeSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties);
	virtual void storeSubSolution(const Step &step, const std::vector<Solution*> &solution, const std::vector<std::pair<ElementType, Property> > &properties);

	virtual void finalize();

	virtual ~ResultStore();

protected:
	ResultStore(const OutputConfiguration &output, const Mesh *mesh, MeshInfo::InfoMode mode = MeshInfo::EMPTY);

	virtual std::vector<std::string> store(const std::string &name, const Step &step, const MeshInfo *meshInfo);
	virtual std::string store(const std::string &name, const RegionData &regionData) =0;
	virtual std::string linkClusters(const std::string &root, const std::string &name, const RegionData &regionData) =0;
	virtual void linkSteps(const std::string &name, const std::vector<std::pair<Step, std::vector<std::string> > > &steps) =0;

	const Mesh *_mesh;

	MeshInfo::InfoMode _mode;
	MeshInfo* _meshInfo;
	std::vector<std::pair<Step, std::vector<std::string> > > _solutions;
	std::vector<std::pair<Step, std::vector<std::string> > > _settings;
	std::vector<std::pair<Step, std::vector<std::string> > > _FETIdata;

private:
	void prepare();
	void storeElementInfo(const Step &step);
	void storeFixPoints(const Step &step);
	void storeCorners(const Step &step);
	void storeDirichlet(const Step &step, const Instance &instance);
	void storeLambdas(const Step &step, const Instance &instance);


};

}

#endif /* SRC_OUTPUT_RESULTSTORE_H_ */
