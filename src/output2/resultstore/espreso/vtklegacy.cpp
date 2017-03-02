
#include <fstream>

#include "../vtklegacy.h"

#include "../../../basis/logging/logging.h"

using namespace espreso::output;

VTKLegacy::VTKLegacy(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path)
{
	preprocessing();
}

void VTKLegacy::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data)
{
	std::ofstream os;

	os.open((name + ".vtk").c_str(), std::ios::out | std::ios::trunc);
	os << "# vtk DataFile Version 4.0\n";
	os << "ESPRESO output\n";
	os << "ASCII\n";
	os << "\n";

	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << coordinates.size() / 3 << " float\n";

	for (size_t i = 0; i < coordinates.size(); i += 3) {
		os << coordinates[i] << " " << coordinates[i + 1] << " " << coordinates[i + 2] << "\n";
	}
	os << "\n";

	os << "CELLS " << elementsTypes.size() << " " << elementsNodes.size() + elements.size()<< "\n";
	for (size_t e = 0, p = 0; e < elementsTypes.size(); p += elementsNodes[e++]) {
		os << elementsNodes[e] << " ";
		for (size_t n = 0; n < elementsNodes[e]; n++) {
			os << elements[p + n] << " ";
		}
		os << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << elementsTypes.size() << "\n";
	for (size_t e = 0, p = 0; e < elementsTypes.size(); e++) {
		os << elementsTypes[e] << "\n";
	}
	os << "\n";

	os << "CELL_DATA " << elementsTypes.size() << "\n";
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		if (it->second->size() % elementsTypes.size() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong integer elements data size.";
		}

		size_t scalarSize = it->second->size() / elementsTypes.size();
		os << "SCALARS " << it->first << " int " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		for (size_t e = 0, i = 0; e < elementsTypes.size(); e++) {
			for (size_t s = 0; s < scalarSize; s++, i++) {
				os << (*it->second)[i] << " ";
			}
			os << "\n";
		}
		os << "\n";
	}

	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		if (it->second->size() % elementsTypes.size() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong double elements data size.";
		}

		size_t scalarSize = it->second->size() / elementsTypes.size();
		os << "SCALARS " << it->first << " int " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		for (size_t e = 0, i = 0; e < elementsTypes.size(); e++) {
			for (size_t s = 0; s < scalarSize; s++, i++) {
				os << (*it->second)[i] << " ";
			}
			os << "\n";
		}
		os << "\n";
	}

	size_t coordinateSize = coordinates.size() / 3;
	os << "POINT_DATA " << coordinateSize << "\n";
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		if (it->second->size() % coordinateSize != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong integer point data size.";
		}

		size_t scalarSize = it->second->size() / coordinateSize;
		os << "SCALARS " << it->first << " int " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		for (size_t e = 0, i = 0; e < coordinateSize; e++) {
			for (size_t s = 0; s < scalarSize; s++, i++) {
				os << (*it->second)[i] << " ";
			}
			os << "\n";
		}
		os << "\n";
	}

	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		if (it->second->size() % coordinateSize != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong double point data size.";
		}

		size_t scalarSize = it->second->size() / coordinateSize;
		os << "SCALARS " << it->first << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		for (size_t e = 0, i = 0; e < coordinateSize; e++) {
			for (size_t s = 0; s < scalarSize; s++, i++) {
				os << (*it->second)[i] << " ";
			}
			os << "\n";
		}
		os << "\n";
	}
}




