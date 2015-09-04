#ifndef FEADAPTOR_HEADER
#define FEADAPTOR_HEADER

#include "../mesh/src/esmesh.h"
#include "../mesh/src/elements/elements.h"
#include "../mesh/src/elements/element.h"
#include "../mesh/src/structures/mesh.h"

class Attributes;
class Grid;

namespace Adaptor
{
  void Initialize(int numScripts, char* scripts[]);

  void Finalize();

  void CoProcess(	mesh::Mesh *mesh,
									std::vector<std::vector<eslocal> > &l2g_vec, 
									std::vector<double>& grid_points, 
									std::vector<std::vector<double> >& prim_solution, 
									std::vector<float>& decomposition_values, double time,
                  unsigned int timeStep, bool lastTimeStep);
}

#endif
