
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esmesh.h"
#include <string>

namespace esinput {

class Loader {

public:
	virtual void points(std::vector<mesh::Point> &data) =0;

	virtual ~Loader() {};
};

template <class TLoader>
class Input {

public:
	Input(int argc, char** argv): _loader(argc, argv) { };

	void load(mesh::Mesh &mesh)
	{
		_loader.points(mesh._coordinates._points);
	}

private:
	TLoader _loader;
};

}


#endif /* INPUT_LOADER_H_ */
