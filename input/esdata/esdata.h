
#ifndef INPUT_ESDATA_ESDATA_H_
#define INPUT_ESDATA_ESDATA_H_

#include "../loader.h"

namespace espreso {
namespace input {

class Esdata: public ExternalLoader {

public:
	Esdata(const Options &options, int rank, int size);

	void points(Coordinates &coordinates);
	void elements(std::vector<Element*> &elements);
	void faces(Faces &faces) {};
	void boundaryConditions(Coordinates &coordinates);
	void clusterBoundaries(Mesh &mesh, Boundaries &boundaries, std::vector<int> &neighbours);

	void open() {};
	void close() {};

private:
	std::string _path;
	int _rank;
	int _size;
};

}
}


#endif /* INPUT_ESDATA_ESDATA_H_ */
