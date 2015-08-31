#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include "../elements/elements.h"
#include "mesh.h"

namespace mesh {

class Mesh;

class Boundaries
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Boundaries &f);

	Boundaries(Mesh &mesh);

	void resize(size_t size)
	{
		_boundaries.resize(size);
	}

	size_t size()
	{
		return _boundaries.size();
	}

	std::set<eslocal>& operator[](size_t position)
	{
		return _boundaries[position];
	}

	void setCorner(size_t index)
	{
		_corners[index] = true;
	}

	bool isCorner(size_t index) const
	{
		return _corners[index];
	}

	Mesh & mesh() {
		return _mesh;
	}


	// prepare for future improvements
	eslocal index(size_t position) {
		return position;
	}

	template<typename T>
	void create_B1_l(	std::vector < SparseIJVMatrix<T> >      	& B1_local,
						std::vector < SparseIJVMatrix<T> >      	& B0_local,
						std::vector < std::vector <eslocal> >      	& l2g_vec,
						std::vector < std::vector <T> >   & lambda_map_sub_clst,
						std::vector < std::vector <T> >   & lambda_map_sub_B1,
						std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
						std::vector < std::vector <double> > 	& B1_l_duplicity,
						const eslocal domains_num,
						mesh::Boundaries & global_boundaries) ;

	template<typename T>
	void create_B1_g(	std::vector < SparseIJVMatrix<T> >         & B1_global,
						const std::vector < SparseCSRMatrix<T> >   & K_mat,
						std::vector < std::vector <eslocal> >   & lambda_map_sub_clst,
						std::vector < std::vector <eslocal> >   & lambda_map_sub_B1,
						std::vector < std::vector <double> >    & B1_duplicity,
						const eslocal MPIrank,
						const eslocal MPIsize,
						const eslocal subDomPerCluster,
						std::vector < eslocal  > & myNeighClusters,
						mesh::Boundaries & local_boundaries) ;



private:
	Mesh &_mesh;

	/** @brief Keeps mapping of nodes to mesh parts. */
	std::vector<std::set<eslocal> > _boundaries;

	/** @brief Keeps information whether a point is the corner. */
	std::vector<bool> _corners;
};

}

#include "boundaries.hpp"

#endif /* BOUNDARIES_H_ */
