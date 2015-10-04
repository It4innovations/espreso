
#ifndef INSTANCE_H_
#define INSTANCE_H_

#include "configuration.h"
#include "esmesh.h"
//#include "espmcube.h"

class Instance {

public:
	Instance(int rank, int size);

	const mesh::Mesh& mesh() const
	{
		return _mesh;
	}

	const mesh::Boundaries& localBoundaries() const
	{
		return _localBoundaries;
	}

	const mesh::Boundaries& globalBoundaries() const
	{
		return _globalBoundaries;
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
	mesh::Boundaries _localBoundaries;
	mesh::Boundaries _globalBoundaries;
};


#endif /* INSTANCE_H_ */
