
#ifndef SRC_OUTPUT_OUTPUT_H_
#define SRC_OUTPUT_OUTPUT_H_

#include "esinfo/stepinfo.h"

#include <string>

namespace espreso {

class DirectOutputExecutor;
class AsyncOutputExecutor;

class OutputWriter {
public:
	virtual bool storeStep() { return true; }

	virtual void updateMesh() =0;
	virtual void updateMonitors(step::TYPE type) =0;

	virtual void updateSolution(const step::Step &step, const step::Time &time) =0;
	virtual void updateSolution(const step::Step &step, const step::Frequency &frequency) =0;

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
	void updateMonitors(step::TYPE type);

	void updateSolution(const step::Step &step, const step::Time &time);
	void updateSolution(const step::Step &step, const step::Frequency &frequency);

	void suppress();
	void permit();

protected:
	DirectOutputExecutor *_direct;
	AsyncOutputExecutor *_async;
};

}

#endif /* SRC_OUTPUT_OUTPUT_H_ */
