
#include "vtk.h"

using namespace esoutput;

void VTK_Surface::coordinatesDisplacement(const std::vector<std::vector<double> > &displacement)
{
	size_t size = 0;
	for (size_t p = 0; p < _surface.parts(); p++) {
		size += _surface.coordinates().localSize(p);
	}

	_vtk << "\n";
	_vtk << "POINT_DATA " << size << "\n";
	_vtk << "SCALARS displacements float 3\n";
	_vtk << "LOOKUP_TABLE default\n";
	for (size_t p = 0; p < displacement.size(); p++) {
		const std::vector<eslocal> &full = _full.coordinates().localToCluster(p);
		const std::vector<eslocal> &surface = _surface.coordinates().localToCluster(p);
		size_t j = 0;
		for (size_t i = 0; i < surface.size(); i++) {
			while (full[j] < surface[i]) {
				j++;
			}
			_vtk << displacement[p][3 * j + 0] << " ";
			_vtk << displacement[p][3 * j + 1] << " ";
			_vtk << displacement[p][3 * j + 2] << "\n";
		}
	}
	_vtk << "\n";
}



