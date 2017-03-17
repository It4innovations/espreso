
#include "../vtklegacy.h"

#include <fstream>

#include "../../regiondata.h"

#include "../../../mesh/structures/elementtypes.h"
#include "../../../assembler/solution.h"
#include "../../../basis/logging/logging.h"

using namespace espreso::output;

VTKLegacy::VTKLegacy(const OutputConfiguration &output, const Mesh *mesh, const std::string &path)
: ResultStore(output, mesh, path)
{

}

static void storeMesh(std::ofstream &os, const espreso::output::RegionData &regionData)
{
	os << "# vtk DataFile Version 4.0\n";
	os << "ESPRESO output\n";
	os << "ASCII\n";
	os << "\n";

	os << "DATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << regionData.coordinates.size() / 3 << " float\n";

	for (size_t i = 0; i < regionData.coordinates.size(); i += 3) {
		os << regionData.coordinates[i] << " " << regionData.coordinates[i + 1] << " " << regionData.coordinates[i + 2] << "\n";
	}
	os << "\n";

	os << "CELLS " << regionData.elementsTypes.size() << " " << regionData.elementsNodes.size() + regionData.elements.size()<< "\n";
	for (size_t e = 0, offset = 0; e < regionData.elementsTypes.size(); offset = regionData.elementsNodes[e++]) {
		os << regionData.elementsNodes[e] - offset << " ";
		for (size_t n = 0; n < regionData.elementsNodes[e] - offset; n++) {
			os << regionData.elements[offset + n] << " ";
		}
		os << "\n";
	}
	os << "\n";

	os << "CELL_TYPES " << regionData.elementsTypes.size() << "\n";
	for (size_t e = 0, p = 0; e < regionData.elementsTypes.size(); e++) {
		os << regionData.elementsTypes[e] << "\n";
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

void VTKLegacy::store(const std::string &name, const RegionData &regionData)
{
	std::ofstream os;
	os.open((name + ".vtk").c_str(), std::ios::out | std::ios::trunc);

	storeMesh(os, regionData);

	os << "CELL_DATA " << regionData.elementsTypes.size() << "\n";
	for (auto it = regionData.data.elementDataInteger.begin(); it != regionData.data.elementDataInteger.end(); ++it) {
		os << "SCALARS " << it->first << " int " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (auto it = regionData.data.elementDataDouble.begin(); it != regionData.data.elementDataDouble.end(); ++it) {
		os << "SCALARS " << it->first << " int " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (size_t i = 0; i < regionData.solutions.size(); i++) {
		if (regionData.solutions[i]->eType != ElementType::ELEMENTS) {
			continue;
		}
		size_t scalarSize = 0;
		for (size_t p = 0; p < regionData.solutions[i]->data.size(); p++) {
			scalarSize += regionData.solutions[i]->data[p].size();
		}
		if (scalarSize % regionData.elementsTypes.size() != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong solution elements data size.";
		}
		scalarSize /= regionData.elementsTypes.size();

		os << "SCALARS " << regionData.solutions[i]->name << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, regionData.solutions[i]->data, scalarSize);
	}

	size_t coordinateSize = regionData.coordinates.size() / 3;
	os << "POINT_DATA " << coordinateSize << "\n";
	for (auto it = regionData.data.pointDataInteger.begin(); it != regionData.data.pointDataInteger.end(); ++it) {
		os << "SCALARS " << it->first << " int " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (auto it = regionData.data.pointDataDouble.begin(); it != regionData.data.pointDataDouble.end(); ++it) {
		os << "SCALARS " << it->first << " double " << it->second.first << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, *it->second.second, it->second.first);
	}

	for (size_t i = 0; i < regionData.solutions.size(); i++) {
		if (regionData.solutions[i]->eType != ElementType::NODES) {
			continue;
		}
		size_t scalarSize = 0;
		for (size_t p = 0; p < regionData.solutions[i]->data.size(); p++) {
			scalarSize += regionData.solutions[i]->data[p].size();
		}
		if (scalarSize % (regionData.coordinates.size() / 3) != 0) {
			ESINFO(ERROR) << "ESPRESO internal error: wrong solution elements data size.";
		}
		scalarSize /= (regionData.coordinates.size() / 3);

		os << "SCALARS " << regionData.solutions[i]->name << " double " << scalarSize << "\n";
		os << "LOOKUP_TABLE default\n";
		storeData(os, regionData.solutions[i]->data, scalarSize);
	}
}




