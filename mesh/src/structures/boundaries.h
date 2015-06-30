#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include "../elements/elements.h"
#include "mesh.h"

namespace mesh {

class Boundaries
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Boundaries &f);

	Boundaries(const Mesh &mesh);

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

	// prepare for future improvements
	eslocal index(size_t position) {
		return position;
	}

	template<typename T>
	void create_B1_l(	std::vector < SparseIJVMatrix >      	& B1_local,
						std::vector < SparseIJVMatrix >      	& B0_local,
						std::vector < std::vector <T> >      	& l2g_vec,
						std::vector < std::vector <eslocal> >   & lambda_map_sub_clst,
						std::vector < std::vector <eslocal> >   & lambda_map_sub_B1,
						std::vector < std::vector <eslocal> >	& lambda_map_sub_B0,
						std::vector < std::vector <double> > 	& B1_l_duplicity,
						const eslocal domains_num) ;

	template<typename T>
	void Boundaries::create_B1_g(	std::vector < SparseIJVMatrix >         & B1_global,
									const std::vector < SparseCSRMatrix >   & K_mat,
									std::vector < std::vector <eslocal> >   & lambda_map_sub_clst,
									std::vector < std::vector <eslocal> >   & lambda_map_sub_B1,
									std::vector < std::vector <double> >    & B1_duplicity,
									const eslocal MPIrank,
									const eslocal MPIsize,
									const eslocal subDomPerCluster,
									std::vector < eslocal  > & myNeighClusters) ;



private:
	const Mesh &_mesh;

	/** @brief Keeps mapping of nodes to mesh parts. */
	std::vector<std::set<eslocal> > _boundaries;
	//std::vector<eslocal> _mapping;
};

}

#include "boundaries.hpp"

#endif /* BOUNDARIES_H_ */
