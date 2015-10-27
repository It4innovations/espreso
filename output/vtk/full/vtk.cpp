
#include "vtk.h"

using namespace esoutput;


void VTK_Full::coordinatesDisplacement(const std::vector<std::vector<double> > &displacement)
{
	size_t size = 0;
	for (size_t p = 0; p < displacement.size(); p++) {
		size += displacement[p].size() / 3;
	}

	_vtk << "\n";
	_vtk << "POINT_DATA " << size << "\n";
	_vtk << "SCALARS displacements float 3\n";
	_vtk << "LOOKUP_TABLE default\n";
	for (size_t p = 0; p < displacement.size(); p++) {
		for (size_t i = 0; i < displacement[p].size() / 3; i++) {
			_vtk << displacement[p][3 * i + 0] << " ";
			_vtk << displacement[p][3 * i + 1] << " ";
			_vtk << displacement[p][3 * i + 2] << "\n";
		}
	}

//	size_t size = 0;
//	for (size_t p = 0; p < displacement.size(); p++) {
//		size += displacement[p].size() / 1;
//	}
//
//	_vtk << "\n";
//	_vtk << "POINT_DATA " << size << "\n";
//	_vtk << "SCALARS displacements float 1\n";
//	_vtk << "LOOKUP_TABLE default\n";
//	for (size_t p = 0; p < displacement.size(); p++) {
//		for (size_t i = 0; i < displacement[p].size(); i++) {
//			_vtk << displacement[p][i] << "\n";
//		}
//	}

	_vtk << "\n";
}



