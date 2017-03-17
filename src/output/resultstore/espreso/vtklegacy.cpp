
#include "../vtklegacy.h"

#include <fstream>

#include "../../../mesh/structures/elementtypes.h"
#include "../../../assembler/solution.h"
#include "../../../basis/logging/logging.h"

using namespace espreso::output;

VTKLegacy::VTKLegacy(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path)
{

}

static void storeMesh(std::ofstream &os, const MeshInfo *regionInfo)
{
	os << "# vtk DataFile Version 4.0\n";
	os << "ESPRESO output\n";
	os << "ASCII\n";
	os << "\n";

	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << regionInfo->coordinates.size() / 3 << " float\n";

	for (size_t i = 0; i < regionInfo->coordinates.size(); i += 3) {
		os << regionInfo->coordinates[i] << " " << regionInfo->coordinates[i + 1] << " " << regionInfo->coordinates[i + 2] << "\n";
	}
	os << "\n";

	os << "CELLS " << regionInfo->elementsTypes.size() << " " << regionInfo->elementsNodes.size() + regionInfo->elements.size()<< "\n";
	for (size_t e = 0, offset = 0; e < regionInfo->elementsTypes.size(); offset = regionInfo->elementsNodes[e++]) {
		os << regionInfo->elementsNodes[e] - offset << " ";
		for (size_t n = 0; n < regionInfo->elementsNodes[e] - offset; n++) {
			os << regionInfo->elements[offset + n] << " ";
		}
		os << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << regionInfo->elementsTypes.size() << "\n";
	for (size_t e = 0, p = 0; e < regionInfo->elementsTypes.size(); e++) {
		os << regionInfo->elementsTypes[e] << "\n";
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

void VTKLegacy::store(const std::string &name, const MeshInfo *regionInfo)
{
	std::ofstream os;
	os.open((name + ".vtk").c_str(), std::ios::out | std::ios::trunc);

	storeMesh(os, regionInfo);

	os << "CELL_DATA " << regionInfo->elementsTypes.size() << "\n";
	for (auto it = regionInfo->data.elementDataInteger.begin(); it != regionInfo->data.elementDataInteger.end(); ++it) {
		os << "SCALARS " << it->first << " int " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (auto it = regionInfo->data.elementDataDouble.begin(); it != regionInfo->data.elementDataDouble.end(); ++it) {
		os << "SCALARS " << it->first << " int " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (size_t i = 0; i < regionInfo->solutions.size(); i++) {
		if (regionInfo->solutions[i]->eType != ElementType::ELEMENTS) {
			continue;
		}
		size_t scalarSize = 0;
		for (size_t p = 0; p < regionInfo->solutions[i]->data.size(); p++) {
			scalarSize += regionInfo->solutions[i]->data[p].size();
		}
		if (scalarSize % regionInfo->elementsTypes.size() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong solution elements data size.";
		}
		scalarSize /= regionInfo->elementsTypes.size();

		os << "SCALARS " << regionInfo->solutions[i]->name << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, regionInfo->solutions[i]->data, scalarSize);
	}

	size_t coordinateSize = regionInfo->coordinates.size() / 3;
	os << "POINT_DATA " << coordinateSize << "\n";
	for (auto it = regionInfo->data.pointDataInteger.begin(); it != regionInfo->data.pointDataInteger.end(); ++it) {
		os << "SCALARS " << it->first << " int " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (auto it = regionInfo->data.pointDataDouble.begin(); it != regionInfo->data.pointDataDouble.end(); ++it) {
		os << "SCALARS " << it->first << " double " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (size_t i = 0; i < regionInfo->solutions.size(); i++) {
		if (regionInfo->solutions[i]->eType != ElementType::NODES) {
			continue;
		}
		size_t scalarSize = 0;
		for (size_t p = 0; p < regionInfo->solutions[i]->data.size(); p++) {
			scalarSize += regionInfo->solutions[i]->data[p].size();
		}
		if (scalarSize % (regionInfo->coordinates.size() / 3) != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong solution elements data size.";
		}
		scalarSize /= (regionInfo->coordinates.size() / 3);

		os << "SCALARS " << regionInfo->solutions[i]->name << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, regionInfo->solutions[i]->data, scalarSize);
	}
}




