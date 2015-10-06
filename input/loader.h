
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esmesh.h"
#include <string>

namespace esinput {

class ExternalLoader {

public:
	void load(mesh::Mesh &mesh)
	{
		points(mesh._coordinates);
		elements(mesh._elements);
	}

	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements) = 0;

	virtual ~ExternalLoader() {};
};

class InternalLoader {

public:
	void load(mesh::Mesh &mesh)
	{
		points(mesh._coordinates);
		elements(mesh._elements, mesh._partPtrs);
		for (size_t i = 0; i < mesh._partPtrs.size() - 1; i++) {
			mesh.computeLocalIndices(i);
		}
	}

	virtual void points(mesh::Coordinates &coordinates) = 0;
	virtual void elements(std::vector<mesh::Element*> &elements, std::vector<eslocal> &parts) = 0;

	virtual ~InternalLoader() {};
};

template <class TLoader>
class Loader {

public:
	Loader(int argc, char** argv, int rank, int size): _loader(argc, argv, rank, size) { };

	void load(mesh::Mesh &mesh)
	{
		_loader.load(mesh);
	}

private:
	TLoader _loader;
};

}


#endif /* INPUT_LOADER_H_ */
