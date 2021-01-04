
#ifndef SRC_OUTPUT_OUTPUT_H_
#define SRC_OUTPUT_OUTPUT_H_

#include <string>
#include <vector>

namespace espreso {

class OutputExecutor;

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
	std::vector<OutputExecutor*> _executors;
};

}

#endif /* SRC_OUTPUT_OUTPUT_H_ */
