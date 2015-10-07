
#ifndef INPUT_ESDATA_ESDATA_H_
#define INPUT_ESDATA_ESDATA_H_

#include "../loader.h"

namespace esinput {

class Esdata: public ExternalLoader {

public:
	Esdata(int argc, char** argv, int rank, int size);

	void points(mesh::Coordinates &coordinates);
	void elements(std::vector<mesh::Element*> &elements);

private:
	std::string _path;
	int _rank;
	int _size;
};

}


#endif /* INPUT_ESDATA_ESDATA_H_ */
