
#ifndef ASSEMBLER_GLUING_GLUING_H_
#define ASSEMBLER_GLUING_GLUING_H_

#include "../assembler.h"

namespace assembler {

// rename to equality constrains
template <class TInput>
class Gluing: public Assembler<TInput> {

public:
	virtual ~Gluing() {};

protected:
	Gluing(TInput &input);

	// computeDirichlet
	// computeGluing
	// computeRigidBodyModes

	void computeSubdomainGluing();
	void computeClusterGluing(std::vector<size_t> &rows);

	void local_B1_global_resize();
	void get_myBorderDOFs_from_mesh();

	virtual size_t DOFs() = 0;

	// Matrices for Linear Solver
	// B0 B1
	std::vector<SparseMatrix> _localB, _globalB;

	std::vector<eslocal> 			   _neighClusters;

	// Description ??
	std::vector<std::vector<eslocal> > _lambda_map_sub_B0;

	std::vector<std::vector<eslocal> > _lambda_map_sub_B1;
	std::vector<std::vector<eslocal> > _lambda_map_sub_clst;
	std::vector<std::vector<double> >  _B1_duplicity;
	std::vector<std::vector<double> >  _vec_c;

private:
	// used for computation
	std::vector<SparseIJVMatrix<eslocal> > _B0;				// local matrix for HFETI corners

	std::vector<SparseIJVMatrix<eslocal> > _B1;				// final B1 - will be created by appending the _B1_dir + _B1_gl_loc + _B1_gl_glob + _B1_gl_glob_red

	std::vector<SparseIJVMatrix<eslocal> > _B1_dir; 		// local B1 for dirichlet
	std::vector<SparseIJVMatrix<eslocal> > _B1_gl_loc;		// local B1 for gluing
	std::vector<SparseIJVMatrix<eslocal> > _B1_gl_glob;		// global B1 for gluing
	std::vector<SparseIJVMatrix<eslocal> > _B1_gl_glob_red; // global B1 for gluing - redundant lagrange multipliers

	// Per-cluster variables
	std::vector < esglobal > _myBorderDOFs;  				// DOFs on the border of the cluster but not domains
	std::vector < esglobal > _myBorderDOFs_sp;  			// DOFs on the border of the cluster and domains

	esglobal 				 total_number_of_B1_l_rows;


};

}

#include "gluing.hpp"


#endif /* ASSEMBLER_GLUING_GLUING_H_ */
