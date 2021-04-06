

#ifndef SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_H_
#define SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_H_

#include "visualization.h"
#include "writer/ensightwriter.h"
#include <string>
#include <functional>

namespace espreso {

class ElementsRegionStore;
class ElementsInterval;

class EnSightGold: public Visualization {

	class FTT;
public:
	EnSightGold(const Mesh &mesh, bool withDecomposition = true);
	~EnSightGold();

	void updateMesh();
	void updateSolution();

protected:
	std::string dataname(const NamedData *data, int d);
	void casefile();
	void geometry();
	void ndata(const NamedData *data);
	void edata(const NamedData *data);
	void decomposition();

	FTT *_ftt;

	bool _withDecomposition;
	std::string _geometry;
	std::string _fixedDataPath;
	std::vector<std::string> _variables;
	std::vector<double> _times;
	EnsightBinaryWriter _writer;
};

class EnSightGold::FTT: public EnSightGold {

public:
	FTT(EnSightGold *parent);
};

}

#endif /* SRC_OUTPUT_VISUALIZATION_ENSIGHTGOLD_H_ */
