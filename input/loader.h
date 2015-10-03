
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esmesh.h"
#include <string>

namespace esinput {

template <class TLoader>
class Loader {

public:
	Loader(const std::string &configuration): _loader(configuration) { };

	void fillMesh(mesh::Mesh &mesh)
	{
		mesh._elements.size();
	}

	virtual ~Loader() {};

private:
	TLoader _loader;
};

}





#endif /* INPUT_LOADER_H_ */
