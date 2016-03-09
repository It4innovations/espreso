
#include "vtk.h"
#include "esconfig.h"

using namespace esoutput;

void VTK_Full::coordinatesDisplacement(const std::vector<std::vector<double> > &displacement, size_t dofs)
{

	size_t size = 0;
	for (size_t p = 0; p < _mesh.parts(); p++) {
		size += _mesh.coordinates().localSize(p);
	}

	_vtk << "\n";
	_vtk << "POINT_DATA " << size << "\n";
	_vtk << "SCALARS displacements float " << dofs << "\n";
	_vtk << "LOOKUP_TABLE default\n";
	for (size_t p = 0; p < displacement.size(); p++) {
		for (size_t i = 0; i < displacement[p].size() / dofs; i++) {
			for (size_t d = 0; d < dofs; d++) {
				_vtk << displacement[p][dofs * i + d] << " ";
			}
			_vtk << "\n";
		}
	}
}



