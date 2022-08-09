
#ifndef SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_VOLUME_H_
#define SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_VOLUME_H_

#include "visualization.h"
#include "writer/ensightwriter.h"
#include <string>
#include <functional>

namespace espreso {

class EnSightGoldVolume: public Visualization {

	class FTT;
public:
	EnSightGoldVolume();
	~EnSightGoldVolume();

	void updateMesh();
	void updateSolution();

protected:
	std::string dataname(const NamedData *data, int d);
	void casefile();
	void geometry();
	int ndata(const NamedData *data);
	int edata(const NamedData *data);
	void decomposition();

	int _step;
	std::string _geometry;
	std::string _fixedDataPath;
	std::vector<std::string> _variables;
	std::vector<double> _times;
	EnsightBinaryOutputWriter _writer;
};

}

#endif /* SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_VOLUME_H_ */
