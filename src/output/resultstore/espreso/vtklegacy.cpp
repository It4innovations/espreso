
#include "../vtklegacy.h"

#include <fstream>

#include "../../../mesh/structures/elementtypes.h"
#include "../../../assembler/solution.h"
#include "../../../basis/logging/logging.h"

using namespace espreso::output;

VTKLegacy::VTKLegacy(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path)
{
	preprocessing();
}

static void storeMesh(std::ofstream &os, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements)
{
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
	for (size_t e = 0, offset = 0; e < elementsTypes.size(); offset = elementsNodes[e++]) {
		os << elementsNodes[e] - offset << " ";
		for (size_t n = 0; n < elementsNodes[e] - offset; n++) {
			os << elements[offset + n] << " ";
		}
		os << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << elementsTypes.size() << "\n";
	for (size_t e = 0, p = 0; e < elementsTypes.size(); e++) {
		os << elementsTypes[e] << "\n";
	}
	os << "\n";
}

template <typename Ttype>
static void storeData(std::ofstream &os, const std::vector<Ttype> &data, size_t scalarSize)
{
	for (size_t i = 0; i < data.size(); i++) {
		os << data[i] << " ";
		if (i && i % scalarSize == 0) {
			os << "\n";
		}
	}
	os << "\n";
}

void VTKLegacy::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const DataArrays &data)
{
	std::ofstream os;
	os.open((name + ".vtk").c_str(), std::ios::out | std::ios::trunc);

	storeMesh(os, coordinates, elementsTypes, elementsNodes, elements);

	os << "CELL_DATA " << elementsTypes.size() << "\n";
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		if (it->second->size() % elementsTypes.size() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong integer elements data size of " << it->first;
		}

		size_t scalarSize = it->second->size() / elementsTypes.size();
		os << "SCALARS " << it->first << " int " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second, scalarSize);
	}

	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		if (it->second->size() % elementsTypes.size() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong double elements data size of " << it->first;
		}

		size_t scalarSize = it->second->size() / elementsTypes.size();
		os << "SCALARS " << it->first << " int " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second, scalarSize);
	}

	size_t coordinateSize = coordinates.size() / 3;
	os << "POINT_DATA " << coordinateSize << "\n";
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		if (it->second->size() % coordinateSize != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong integer point data size of " << it->first;
		}

		size_t scalarSize = it->second->size() / coordinateSize;
		os << "SCALARS " << it->first << " int " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second, scalarSize);
	}

	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		if (it->second->size() % coordinateSize != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong double point data size of " << it->first;
		}

		size_t scalarSize = it->second->size() / coordinateSize;
		os << "SCALARS " << it->first << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second, scalarSize);
	}
}

template <typename Ttype>
static void storeData(std::ofstream &os, const std::vector<std::vector<Ttype> > &data, size_t scalarSize)
{
	for (size_t p = 0; p < data.size(); p++) {
		for (size_t i = 0; i < data[p].size(); i++) {
			os << data[p][i] << " ";
			if ((i || p) && i % scalarSize == 0) {
				os << "\n";
			}
		}
	}
	os << "\n";
}

void VTKLegacy::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution)
{
	std::ofstream os;
	os.open((name + ".vtk").c_str(), std::ios::out | std::ios::trunc);

	storeMesh(os, coordinates, elementsTypes, elementsNodes, elements);

	os << "CELL_DATA " << elementsTypes.size() << "\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType != ElementType::ELEMENTS) {
			continue;
		}
		size_t scalarSize = 0;
		for (size_t p = 0; p < solution[i]->data.size(); p++) {
			scalarSize += solution[i]->data[p].size();
		}
		if (scalarSize % elementsTypes.size() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong solution elements data size.";
		}
		scalarSize /= elementsTypes.size();

		os << "SCALARS " << solution[i]->name << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, solution[i]->data, scalarSize);
	}

	os << "POINT_DATA " << coordinates.size() / 3 << "\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType != ElementType::NODES) {
			continue;
		}
		size_t scalarSize = 0;
		for (size_t p = 0; p < solution[i]->data.size(); p++) {
			scalarSize += solution[i]->data[p].size();
		}
		if (scalarSize % (coordinates.size() / 3) != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong solution elements data size.";
		}
		scalarSize /= (coordinates.size() / 3);

		os << "SCALARS " << solution[i]->name << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, solution[i]->data, scalarSize);
	}
}




