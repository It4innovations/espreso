#ifndef ADAPTOR_H
#define ADAPTOR_H

#include <iostream>
#include <vector>

#include "mpi.h"

#include "esmesh.h"

#define MAX_NUMBER_OF_NODES 20

class Attributes;
class Grid;

namespace Adaptor
{
  void Initialize(int numScripts, char* scripts[]);

  void Finalize();

  void CoProcess(
			mesh::Mesh *mesh,
			std::vector< std::vector<eslocal> > &l2g_vec,
			std::vector< std::vector<double> > & prim_solution,
			double time,
			unsigned int timeStep,
			bool lastTimeStep);
}

#endif
