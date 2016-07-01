
#include "vtk.h"

using namespace espreso::output;

void VTK_Surface::coordinatesDisplacement(const std::vector<std::vector<double> > &displacement)
{
	size_t DOFs = displacement[0].size() / _surface.coordinates().localSize(0);
	size_t size = 0;
	for (size_t p = 0; p < _surface.parts(); p++) {
		size += _surface.coordinates().localSize(p);
	}

	_vtk << "\n";
	_vtk << "POINT_DATA " << size << "\n";
	_vtk << "SCALARS displacements float " << DOFs << "\n";
	_vtk << "LOOKUP_TABLE default\n";
	for (size_t p = 0; p < displacement.size(); p++) {
		const std::vector<eslocal> &full = _full.coordinates().localToCluster(p);
		const std::vector<eslocal> &surface = _surface.coordinates().localToCluster(p);
		size_t j = 0;
		for (size_t i = 0; i < surface.size(); i++) {
			while (full[j] < surface[i]) {
				j++;
			}
			for (size_t d = 0; d < DOFs; d++) {
				_vtk << displacement[p][DOFs * j + d] << " ";
			}
			_vtk << "\n";
		}
	}
	_vtk << "\n";
}



