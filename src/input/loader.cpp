
#include "loader.h"

#include "../mesh/settings/property.h"
#include "../config/description.h"

#include "ansys/ansys.h"
#include "openfoam/openfoam.h"
#include "esdata/esdata.h"
#include "generator/generator.h"

using namespace espreso::input;


void Loader::load(const GlobalConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	switch (configuration.input) {
	case INPUT::WORKBENCH:
		AnsysWorkbench::load(configuration, mesh, index, size);
		break;
	case INPUT::OPENFOAM:
		OpenFOAM::load(configuration, mesh, index, size);
		break;
	case INPUT::ESDATA:
		Esdata::load(configuration, mesh, index, size);
		break;
	case INPUT::GENERATOR:
		Generator::generate(configuration, mesh, index, size);
		break;
	}
}

void Loader::fill()
{
	auto intervalStats = [] (const std::vector<eslocal> &data) {
		std::vector<size_t> sizes(data.size() - 1);
		for (size_t p = 0; p < data.size() - 1; p++) {
			sizes[p] = data[p + 1] - data[p];
		}
		return Info::averageValues(sizes);
	};

	TimeEval measurement("Mesh loader"); measurement.totalTime.startWithBarrier();

	open();

	// LOAD POINTS
	TimeEvent tPoints("coordinates"); tPoints.start();
	points(mesh._coordinates);
	std::sort(mesh._coordinates._globalMapping.begin(), mesh._coordinates._globalMapping.end());
	tPoints.end(); measurement.addEvent(tPoints);
	ESINFO(OVERVIEW) << "Coordinates loaded - total number of nodes: " << Info::sumValue(mesh.coordinates().clusterSize());

	// LOAD ELEMENTS
	TimeEvent tElements("elements"); tElements.start();
	elements(mesh._elements, mesh._faces, mesh._edges);
	materials(mesh._materials);

	mesh.fillNodesFromCoordinates();
	mesh.fillParentElementsToNodes();
	tElements.end(); measurement.addEvent(tElements);
	ESINFO(OVERVIEW) << "Elements loaded - total number of elements: " << Info::sumValue(mesh.elements().size());

	// LOAD NEIGHBOURS
	TimeEvent tClusterBoundaries("cluster boundaries"); tClusterBoundaries.start();
	neighbours(mesh._nodes, mesh._neighbours, mesh._faces, mesh._edges);
	tClusterBoundaries.end(); measurement.addEvent(tClusterBoundaries);
	ESINFO(OVERVIEW) << "Neighbours loaded - number of neighbours for each cluster is " << Info::averageValue(mesh.neighbours().size());


	// LOAD BOUNDARY CONDITIONS
	TimeEvent tSettings("settings"); tSettings.start();
	regions(mesh._evaluators, mesh._regions, mesh._elements, mesh._faces, mesh._edges, mesh._nodes);

	for (size_t r = 0; r < mesh._regions.size(); r++) {
		ESINFO(OVERVIEW) << "Loaded region '" << mesh._regions[r].name << "' of size " << Info::sumValue(mesh._regions[r].elements.size());
	}

	mesh.fillFacesParents();
	mesh.fillEdgesParents();

	tSettings.end(); measurement.addEvent(tSettings);

	close();

	TimeEvent tPartition("partition"); tPartition.start();
	if (partitiate(mesh._nodes, mesh._partPtrs, mesh._fixPoints, mesh._corners)) { // manual partition -> map elements to the domains
		mesh.mapElementsToDomains();
		mesh.mapFacesToDomains();
		mesh.mapEdgesToDomains();
		mesh.mapNodesToDomains();
		mesh.mapCoordinatesToDomains();
	}

	tPartition.end(); measurement.addEvent(tPartition);
	ESINFO(OVERVIEW) << "Mesh partitioned into " << environment->MPIsize << " * " << mesh.parts() << " = " << mesh.parts() * environment->MPIsize
			<< " parts. There is " << intervalStats(mesh._partPtrs) << " elements in subdomain.";

	measurement.totalTime.endWithBarrier(); measurement.printStatsMPI();
}



