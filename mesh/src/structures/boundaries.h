#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include "../elements/elements.h"
#include "mesh.h"

class Boundaries
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Boundaries &f);

	Boundaries(const Mesh &mesh, const Coordinates &coordinates);
	~Boundaries() { };

	void create_B1_l(	std::vector < SparseIJVMatrix >      & B1_local, 
									std::vector < SparseIJVMatrix >      & B0_local,
									std::vector < std::vector <int> >    & l2g_vec,
									std::vector < std::vector <int> >	 & lambda_map_sub_clst,
									std::vector < std::vector <int> >    & lambda_map_sub_B1, 
									std::vector < std::vector <int> >    & lambda_map_sub_B0, 
									std::vector < std::vector <double> > & B1_l_duplicity,
									const int domains_num, 
									const Mesh &mesh) ;


private:
	/** @brief Reference to coordinates. */
	const Coordinates &_coordinates;

	/** @brief Keeps mapping of nodes to mesh parts. */
	BoundaryNodes _boundaries;
};



#endif /* BOUNDARIES_H_ */
