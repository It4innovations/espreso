
#ifndef SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_
#define SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_

#include "output/output.h"
#include "mesh/store/statisticsstore.h"

#include <utility>
#include <vector>
#include <cstdio>

namespace espreso {

struct ElementsRegionStore;
struct BoundaryRegionStore;
struct NodeData;
struct ElementData;

struct Monitor {
	std::string name;
	std::string property;
	std::string stats;
	int printSize;
	double *data;

	Monitor(): name("---"), property("-"), printSize(5), data(NULL) {}
};

class Monitoring: public OutputWriter {

public:
	bool storeStep(const step::Step &step);

	void updateMesh() {}
	void updateMonitors(const step::Step &step);
	void updateSolution(const step::Step &step, const step::Time &time);
	void updateSolution(const step::Step &step, const step::Frequency &frequency);

	Monitoring();
	~Monitoring();
	static char delimiter;

protected:
	void updateSolution(const step::Step &step);
	void storeSolution(const step::Step &step);

	FILE *_runFile, *_fttFile;


	std::vector<Monitor> _monitors;

	std::vector<std::pair<NodeData*, const ElementsRegionStore*> > _nedata;
	std::vector<std::pair<NodeData*, const BoundaryRegionStore*> > _nbdata;
	std::vector<std::pair<ElementData*, const ElementsRegionStore*> > _edata;

	std::vector<Statistics> _statistics;
};

}



#endif /* SRC_OUTPUT_RESULT_MONITORS_MONITORING_H_ */
