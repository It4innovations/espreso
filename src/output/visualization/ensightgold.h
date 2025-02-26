

#ifndef SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_H_
#define SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_H_

#include "visualization.h"
#include "writer/ensightwriter.h"
#include <string>
#include <functional>

namespace espreso {

struct ElementsRegionStore;
struct ElementsInterval;

class EnSightGold: public Visualization {

	class FTT;
public:
	EnSightGold();
	~EnSightGold();

	void updateMesh();
	void updateMonitors(const step::Step &step);
	void updateSolution(const step::Step &step, const step::Time &time);
	void updateSolution(const step::Step &step, const step::Frequency &frequency);

protected:
	void updateSolution(const step::Step &step);
	std::string dataname(const NamedData *data, int d);
	void casefile();
	void geometry();
	int ndata(const NamedData *data, const step::Step &step);
	int edata(const NamedData *data, const step::Step &step);
	void decomposition();

	FTT *_ftt;

	bool _withIDs;
	int _step;
	std::string _geometry;
	std::string _fixedDataPath;
	std::vector<std::string> _variables;
	std::vector<double> _times;
	EnsightBinaryOutputWriter _writer;
};

class EnSightGold::FTT: public EnSightGold {

public:
	FTT(EnSightGold *parent, const step::Frequency &frequency);
};

}

#endif /* SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_H_ */
