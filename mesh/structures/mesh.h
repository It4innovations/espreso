#ifndef MESH_H_
#define MESH_H_

#include <cstring>
#include <algorithm>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "mkl_spblas.h"
#include "mkl_blas.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"

#include "cilk/cilk.h"

#include "metis.h"

#include "../elements/elements.h"
#include "coordinates.h"
#include "boundaries.h"

#include "esbasis.h"
#include "esconfig.h"

namespace esinput {
class InternalLoader;
class ExternalLoader;
}

namespace mesh {


class Boundaries;

class Mesh
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Mesh &m);
	friend class esinput::InternalLoader;
	friend class esinput::ExternalLoader;


	Mesh();

	const Coordinates& coordinates() const
	{
		return _coordinates;
	}

	Coordinates& coordinates()
	{
		return _coordinates;
	}

	const Boundaries& subdomainBoundaries() const
	{
		return _subdomainBoundaries;
	}

	const Boundaries& clusterBoundaries() const
	{
		return _clusterBoundaries;
	}

	void saveNodeArray(eslocal *nodeArray, size_t part) const;

	void getSurface(Mesh &surface) const;

	~Mesh();

	void partitiate(size_t parts);
	void computeFixPoints(size_t fixPoints);
	void computeCorners(eslocal number, bool vertex, bool edges, bool faces, bool averaging);

	const std::vector<Element*>& getElements() const
	{
		return _elements;
	};

	size_t parts() const
	{
		return _partPtrs.size() - 1;
	}

	const std::vector<eslocal>& getPartition() const
	{
		return _partPtrs;
	}

	const std::vector<std::vector<eslocal> >& getFixPoints() const
	{
		return _fixPoints;
	}

	eslocal getPartNodesCount(eslocal part) const
	{
		return _coordinates.localSize(part);
	}

protected:
	eslocal* getPartition(eslocal first, eslocal last, eslocal parts) const;
	eslocal getCentralNode(eslocal first, eslocal last, eslocal *ePartition, eslocal part, eslocal subpart) const;

	void computeBoundaries();
	void remapElementsToSubdomain();
	void remapElementsToCluster();

	/** @brief Reference to coordinates. */
	Coordinates _coordinates;

	/** @brief Array that stores all elements of the mesh. */
	std::vector<mesh::Element*> _elements;

	/** @brief Elements in part 'i' are from _partPtrs[i] to _partPtrs[i + 1]. */
	std::vector<eslocal> _partPtrs;

	/** @brief Fix points for all parts. */
	std::vector<std::vector<eslocal> > _fixPoints;

	/** @brief Map of points to sub-domains. */
	Boundaries _subdomainBoundaries;

	/** @brief Map of points to clusters. */
	Boundaries _clusterBoundaries;
};

}


#endif /* MESH_H_ */
