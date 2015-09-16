#ifndef ADAPTOR_H
#define ADAPTOR_H

#include <iostream>
#include <vector>

#include "mpi.h"

#include "../mesh/src/esmesh.h"
#include "../mesh/src/elements/elements.h"
#include "../mesh/src/elements/element.h"
#include "../mesh/src/structures/mesh.h"

#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

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
