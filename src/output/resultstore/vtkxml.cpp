
#include "vtkxml.h"

#include <fstream>

#include "../../configuration/environment.h"
#include "../../configuration/output.h"
#include "../../assembler/solution.h"
#include "../../mesh/structures/elementtypes.h"

using namespace espreso::output;

void VTKXML::store(const std::string &name, const RegionInfo *regionInfo)
{
	initWriter(name, regionInfo->coordinates.size() / 3, regionInfo->elementsTypes.size());
	addMesh(regionInfo);
	addData(regionInfo->data, regionInfo->solutions);
	finalizeWriter();
}

void VTKXML::linkClusters(const std::string &root, const std::string &name, const RegionInfo *regionInfo)
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
	for (auto it = regionInfo->data.pointDataInteger.begin(); it != regionInfo->data.pointDataInteger.end(); ++it) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (auto it = regionInfo->data.pointDataDouble.begin(); it != regionInfo->data.pointDataDouble.end(); ++it) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (size_t i = 0; i < regionInfo->solutions.size(); i++) {
		if (regionInfo->solutions[i]->eType == ElementType::NODES) {
			size_t size = 0;
			for (size_t p = 0; p < regionInfo->solutions[i]->data.size(); p++) {
				size += regionInfo->solutions[i]->data[p].size();
			}
			os << "      <PDataArray type=\"Float64\" Name=\"" << regionInfo->solutions[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << regionInfo->solutions[i]->properties << "\"/>\n";
		}
	}
	os << "    </PPointData>\n";
	os << "\n";
	os << "    <PCellData>\n";
	for (auto it = regionInfo->data.elementDataInteger.begin(); it != regionInfo->data.elementDataInteger.end(); ++it) {
		os << "      <PDataArray type=\"Int" << 8 * sizeof(eslocal) << "\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (auto it = regionInfo->data.elementDataDouble.begin(); it != regionInfo->data.elementDataDouble.end(); ++it) {
		os << "      <PDataArray type=\"Float64\" Name=\"" << it->first << "\" format=\"" << format() << "\" NumberOfComponents=\"" << it->second.first << "\"/>\n";
	}
	for (size_t i = 0; i < regionInfo->solutions.size(); i++) {
		if (regionInfo->solutions[i]->eType == ElementType::ELEMENTS) {
			size_t size = 0;
			for (size_t p = 0; p < regionInfo->solutions[i]->data.size(); p++) {
				size += regionInfo->solutions[i]->data[p].size();
			}
			os << "      <PDataArray type=\"Float64\" Name=\"" << regionInfo->solutions[i]->name << "\" format=\"" << format() << "\" NumberOfComponents=\"" << regionInfo->solutions[i]->properties << "\"/>\n";
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



