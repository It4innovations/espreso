
#ifndef INSTANCE_H_
#define INSTANCE_H_

#include "configuration.h"
#include "esmesh.h"
//#include "espmcube.h"

class Instance {

public:
	Instance(int argc, char** argv, int rank, int size);

	const mesh::Mesh& mesh() const
	{
		return _mesh;
	}

	const mesh::Boundaries& localBoundaries() const
	{
		return _mesh.subdomainBoundaries();
	}

	const mesh::Boundaries& globalBoundaries() const
	{
		return _mesh.clusterBoundaries();
	}

	int rank() const
	{
		return _rank;
	}

	int size() const
	{
		return _size;
	}

private:
	int _rank;
	int _size;

	mesh::Mesh _mesh;
};


#endif /* INSTANCE_H_ */
