#ifndef FEADAPTOR_HEADER
#define FEADAPTOR_HEADER


class Attributes;
class Grid;

namespace Adaptor
{
  void Initialize(int numScripts, char* scripts[]);

  void Finalize();

  void CoProcess(std::vector<double>& grid_points, std::vector<unsigned int>& cell_points, std::vector<std::vector<double> >& prim_solution, std::vector<float>& decomposition_values, double time,
                 unsigned int timeStep, bool lastTimeStep);
}

#endif