
#ifndef SRC_OUTPUT_RESULT_RESULTSTORE_H_
#define SRC_OUTPUT_RESULT_RESULTSTORE_H_

#include <string>

namespace espreso {

class Mesh;
class ResultStoreExecutor;
class Dispatcher;

class ResultStoreBase {

public:
	virtual void updateMesh() =0;
	virtual void updateMonitors() =0;
	virtual void updateSolution() =0;

	virtual const Mesh& mesh() const { return _mesh; }

	virtual ~ResultStoreBase() {};

protected:
	ResultStoreBase(const Mesh &mesh);

	void createOutputDirectory();

	const Mesh &_mesh;
	std::string _path;
	std::string _directory;
	std::string _name;
	bool _measure;
};

class ResultStore {

public:
	static void createAsynchronizedStore();
	static void destroyAsynchronizedStore();
	static bool isStoreNode();
	static bool isComputeNode();

	void suppress();
	void permit();

	bool storeStep();

	void updateMesh();
	void updateMonitors();
	void updateSolution();

	ResultStore();
	~ResultStore();

protected:
	bool _allowed;

	ResultStoreExecutor *_async;
	ResultStoreExecutor *_direct;

	static Dispatcher *_dispatcher;
};

}

#endif /* SRC_OUTPUT_RESULT_RESULTSTORE_H_ */
