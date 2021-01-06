
#ifndef SRC_OUTPUT_OUTPUT_H_
#define SRC_OUTPUT_OUTPUT_H_

#include <string>

namespace espreso {

class DirectOutputExecutor;
class AsyncOutputExecutor;

class OutputWriter {
public:
	virtual bool storeStep() { return true; }

	virtual void updateMesh() = 0;
	virtual void updateMonitors() = 0;
	virtual void updateSolution() = 0;

	virtual ~OutputWriter() {};

protected:
	OutputWriter();
	void createOutputDirectory();

	std::string _path, _directory, _name;
	bool _measure, _allowed;
};

class Output: public OutputWriter {
public:
	Output();
	~Output();

	void updateMesh();
	void updateMonitors();
	void updateSolution();

	void suppress();
	void permit();

protected:
	DirectOutputExecutor *_direct;
	AsyncOutputExecutor *_async;
};

}

#endif /* SRC_OUTPUT_OUTPUT_H_ */
