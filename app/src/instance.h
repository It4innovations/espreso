
#ifndef INSTANCE_H_
#define INSTANCE_H_

#include "configuration.h"
#include "esmesh.h"
#include "espmcube.h"

class Instance {

public:
	Instance(const Configuration &configuration, int rank, int size);

	const Configuration& configuration() const
	{
		return _configuration;
	}

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

	Configuration _configuration;

	int _rank;
	int _size;

	mesh::Mesh _mesh;
	mesh::Boundaries _localBoundaries;
	mesh::Boundaries _globalBoundaries;
};


#endif /* INSTANCE_H_ */
