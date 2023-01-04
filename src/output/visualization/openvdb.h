
#ifndef SRC_OUTPUT_VISUALIZATION_OPENVDB_H_
#define SRC_OUTPUT_VISUALIZATION_OPENVDB_H_

#include "visualization.h"
#include "basis/containers/allocators.h"
#include "mesh/store/elementstore.h"
#include "wrappers/mpi/communication.h"
#include "wrappers/openvdb/w.openvdb.h"
#include "wrappers/pthread/w.pthread.h"

#include <string>
#include <queue>

struct SharedVolume;

namespace espreso {

class OpenVDB: public Visualization {
public:
	OpenVDB();
	~OpenVDB();

	void lock();
	void updateMesh();
	void updateMonitors();
	void updateSolution();

protected:
	std::string _filename;
	int _step;
	std::vector<ElementData*> activeVariables;

	struct OpenVDBData {
		SharedVolume *volume;
		int root;
		MPI_Request req;

		OpenVDBData(SharedVolume *volume);
		bool test();
	};

	std::queue<OpenVDBData> _postponed;
	std::vector<int> nranks;
};

}

#endif
