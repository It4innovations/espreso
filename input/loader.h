
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

	virtual void points(mesh::Coordinates &data) = 0;
	virtual void elements(std::vector<mesh::Element*> &data) = 0;

	virtual ~ExternalLoader() {};
};

class InternalLoader {

public:
	void load(mesh::Mesh &mesh)
	{
		points(mesh._coordinates);
		elements(mesh._elements);
	}

	virtual void points(mesh::Coordinates &data) = 0;
	virtual void elements(std::vector<mesh::Element*> &data) = 0;

	virtual ~InternalLoader() {};
};

template <class TLoader>
class Input {

public:
	Input(int argc, char** argv, int rank, int size): _loader(argc, argv, rank, size) { };

	void load(mesh::Mesh &mesh)
	{
		_loader.load(mesh);
	}

private:
	TLoader _loader;
};

}


#endif /* INPUT_LOADER_H_ */
