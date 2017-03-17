
#include "vtkxml.h"

#include <fstream>

#include "../../configuration/environment.h"
#include "../../configuration/output.h"
#include "../../assembler/solution.h"
#include "../../mesh/structures/elementtypes.h"
#include "../regiondata.h"

using namespace espreso::output;

void VTKXML::store(const std::string &name, const RegionData &regionData)
{
	initWriter(name, regionData.coordinates.size() / 3, regionData.elementsTypes.size());
	addMesh(regionData);
	addData(regionData.data, regionData.solutions);
	finalizeWriter();
}

void VTKXML::linkClusters(const std::string &root, const std::string &name, const RegionData &regionData)
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
	for (auto it = regionData.data.pointDataInteger.begin(); it != regionData.data.pointDataInteger.end(); ++it) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (auto it = regionData.data.pointDataDouble.begin(); it != regionData.data.pointDataDouble.end(); ++it) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (size_t i = 0; i < regionData.solutions.size(); i++) {
		if (regionData.solutions[i]->eType == ElementType::NODES) {
			size_t size = 0;
			for (size_t p = 0; p < regionData.solutions[i]->data.size(); p++) {
				size += regionData.solutions[i]->data[p].size();
			}
			os << "      <PDataArray type=\"Float64\" Name=\"" << regionData.solutions[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << regionData.solutions[i]->properties << "\"/>\n";
		}
	}
	os << "    </PPointData>\n";
	os << "\n";
	os << "    <PCellData>\n";
	for (auto it = regionData.data.elementDataInteger.begin(); it != regionData.data.elementDataInteger.end(); ++it) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (auto it = regionData.data.elementDataDouble.begin(); it != regionData.data.elementDataDouble.end(); ++it) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (size_t i = 0; i < regionData.solutions.size(); i++) {
		if (regionData.solutions[i]->eType == ElementType::ELEMENTS) {
			size_t size = 0;
			for (size_t p = 0; p < regionData.solutions[i]->data.size(); p++) {
				size += regionData.solutions[i]->data[p].size();
			}
			os << "      <PDataArray type=\"Float64\" Name=\"" << regionData.solutions[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << regionData.solutions[i]->properties << "\"/>\n";
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

void VTKXML::linkSteps(const std::string &name, const std::vector<std::pair<std::string, Step> > &steps)
{
	std::ofstream os;

	std::string ext = _configuration.collected ? ".vtu" : ".pvtu";

	os.open(name + ".pvd", std::ios::out | std::ios::trunc);

	os << "<?xml version=\"1.0\"?>\n";
	os << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	os << "<Collection>\n";
	for (size_t i = 0; i < steps.size(); i++) {
		os << "  <DataSet timestep=\"" << steps[i].second.currentTime << "\" file=\"" << steps[i].first << name << ext << "\"/>\n";
	}
	os << "</Collection>\n";
	os << "</VTKFile>\n";
}



