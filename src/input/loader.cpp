
#include "loader.h"

#include "ansys/ansys.h"
#include "openfoam/openfoam.h"
#include "generator/generator.h"

#include "../mesh/elements/element.h"
#include "../mesh/structures/mesh.h"
#include "../mesh/structures/region.h"
#include "../mesh/settings/evaluator.h"

#include "../basis/logging/timeeval.h"
#include "../basis/logging/logging.hpp"
#include "../basis/utilities/utils.h"
#include "../configuration/globalconfiguration.h"
#include "espreso/espresobinaryformat.h"

using namespace espreso::input;


void Loader::load(const GlobalConfiguration &configuration, Mesh &mesh, size_t index, size_t size)
{
	switch (configuration.input) {
	case INPUT::WORKBENCH:
		AnsysWorkbench::load(configuration.workbench, mesh, index, size);
		break;
	case INPUT::OPENFOAM:
		OpenFOAM::load(configuration.openfoam, mesh, index, size);
		break;
	case INPUT::ESDATA:
		ESPRESOBinaryFormat::load(configuration.esdata, mesh, index, size);
		break;
	case INPUT::GENERATOR:
		Generator::generate(configuration.generator, mesh, index, size);
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
	points(*mesh._coordinates);
	std::sort(mesh._coordinates->_globalMapping.begin(), mesh._coordinates->_globalMapping.end());
	tPoints.end(); measurement.addEvent(tPoints);
	ESINFO(OVERVIEW) << "Coordinates loaded - total number of nodes: " << Info::sumValue(mesh.coordinates().clusterSize());

	// LOAD ELEMENTS
	TimeEvent tElements("elements"); tElements.start();
	elements(mesh._bodies, mesh._elements, mesh._faces, mesh._edges);
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

	if (Test::report(EXPENSIVE)) {
		mesh.checkNeighbours();
	}

	// LOAD BOUNDARY CONDITIONS
	TimeEvent tSettings("settings"); tSettings.start();
	regions(mesh._evaluators, mesh._regions, mesh._elements, mesh._faces, mesh._edges, mesh._nodes);

	for (size_t r = 0; r < mesh._regions.size(); r++) {
		ESINFO(OVERVIEW) << "Loaded region '" << mesh._regions[r]->name << "' of size " << Info::sumValue(mesh._regions[r]->elements().size());
	}

	for (size_t r = 0; r < mesh._regions.size(); r++) {

		size_t threads = environment->OMP_NUM_THREADS;
		std::vector<size_t> distribution = Esutils::getDistribution(threads, mesh._regions[r]->elements().size());

		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				mesh._regions[r]->elements()[n]->regions().push_back(mesh._regions[r]);
			}
		}

		for (size_t s = 0; s < mesh._regions[r]->settings.size(); s++) {
			for (auto it = mesh._regions[r]->settings[s].begin(); it != mesh._regions[r]->settings[s].end(); ++it) {
				ESINFO(OVERVIEW) << it->first << " loaded for LOAD STEP " << s + 1 << " for region '" << mesh._regions[r]->name << "'";
			}
		}
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



