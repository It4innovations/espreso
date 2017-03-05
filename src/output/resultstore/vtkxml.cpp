
#include "vtkxml.h"

#include <fstream>
#include <algorithm>

#include "../../configuration/environment.h"
#include "../../assembler/solution.h"
#include "../../mesh/structures/elementtypes.h"

using namespace espreso::output;

void VTKXML::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, DataArrays &data)
{
	initWriter(name, coordinates.size() / 3, elementsTypes.size());
	addMesh(coordinates, elementsTypes, elementsNodes, elements);
	addData(coordinates.size() / 3, elementsTypes.size(), data);
	finalizeWriter();
}

void VTKXML::store(const std::string &name, std::vector<double> &coordinates, std::vector<eslocal> &elementsTypes, std::vector<eslocal> &elementsNodes, std::vector<eslocal> &elements, const std::vector<Solution*> &solution)
{
	initWriter(name, coordinates.size() / 3, elementsTypes.size());
	addMesh(coordinates, elementsTypes, elementsNodes, elements);
	addData(coordinates.size() / 3, elementsTypes.size(), solution);
	finalizeWriter();
}

void VTKXML::linkClusters(const std::string &root, const std::string &name, const DataArrays &data)
{
	std::ofstream os;

	os.open(root + name + ".pvtu", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
	os << "\n";
	os << "    <PPoints>\n";
	os << "      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << format() << "\"/>\n";
	os << "    </PPoints>\n";
	os << "\n";
	os << "    <PPointData>\n";
	for (auto it = data.pointDataInteger.begin(); it != data.pointDataInteger.end(); ++it) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second << "\"/>\n";
	}
	for (auto it = data.pointDataDouble.begin(); it != data.pointDataDouble.end(); ++it) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second << "\"/>\n";
	}
	os << "    </PPointData>\n";
	os << "\n";
	os << "    <PCellData>\n";
	for (auto it = data.elementDataInteger.begin(); it != data.elementDataInteger.end(); ++it) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second << "\"/>\n";
	}
	for (auto it = data.elementDataDouble.begin(); it != data.elementDataDouble.end(); ++it) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second << "\"/>\n";
	}
	os << "    </PCellData>\n";
	os << "\n";
	for (int r = 0; r < environment->MPIsize; r++) {
		os << "    <Piece Source=\"" << r << "/" << name << ".vtu\"/>\n";
	}
	os << "  </PUnstructuredGrid>\n";
	os << "</VTKFile>\n";
}

void VTKXML::linkClusters(const std::string &root, const std::string &name, const std::vector<Solution*> &solution, size_t points, size_t cells)
{
	std::ofstream os;

	os.open(root + name + ".pvtu", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n";
	os << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
	os << "\n";
	os << "    <PPoints>\n";
	os << "      <PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"" << format() << "\"/>\n";
	os << "    </PPoints>\n";
	os << "\n";
	os << "    <PPointData>\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::NODES) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= points;
			os << "    <DataArray type=\"Float64\" Name=\"" << solution[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << components  << "\">\n";
			os << "      ";
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { store(os, part); });;
			os << "\n";
			os << "    </DataArray>\n";
		}
	}
	os << "    </PPointData>\n";
	os << "\n";
	os << "    <PCellData>\n";
	for (size_t i = 0; i < solution.size(); i++) {
		if (solution[i]->eType == ElementType::ELEMENTS) {
			size_t components = 0;
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { components += part.size(); });
			components /= cells;
			os << "    <DataArray type=\"Float64\" Name=\"" << solution[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << components << "\">\n";
			os << "      ";
			std::for_each(solution[i]->data.begin(), solution[i]->data.end(), [&] (const std::vector<double> &part) { store(os, part); });;
			os << "\n";
			os << "    </DataArray>\n";
		}
	}
	os << "    </PCellData>\n";
	os << "\n";
	for (int r = 0; r < environment->MPIsize; r++) {
		os << "    <Piece Source=\"" << r << "/" << name << ".vtu\"/>\n";
	}
	os << "  </PUnstructuredGrid>\n";
	os << "</VTKFile>\n";
}

void VTKXML::linkSteps(const std::string &root, const std::string &name, const DataArrays &data)
{

}

void VTKXML::linkSteps(const std::string &root, const std::string &name, const std::vector<Solution*> &solution)
{

}


