
#ifndef INPUT_LOADER_H_
#define INPUT_LOADER_H_

#include "esbasis.h"
#include "esmesh.h"
#include <string>

namespace espreso {
namespace input {

class Loader {

public:
	void fill()
	{
		auto sizeStats = [] (const std::vector<std::vector<eslocal> > &data) {
			std::vector<size_t> sizes(data.size());
			for (size_t p = 0; p < data.size(); p++) {
				sizes[p] = data[p].size();
			}
			return Info::averageValues(sizes);
		};
		auto intervalStats = [] (const std::vector<eslocal> &data) {
			std::vector<size_t> sizes(data.size() - 1);
			for (size_t p = 0; p < data.size() - 1; p++) {
				sizes[p] = data[p + 1] - data[p];
			}
			return Info::averageValues(sizes);
		};

		TimeEval measurement("Mesh loader"); measurement.totalTime.startWithBarrier();

		open();

		TimeEvent tPoints("coordinates"); tPoints.start();
		points(mesh._coordinates);
		tPoints.end(); measurement.addEvent(tPoints);
		ESINFO(OVERVIEW) << "Coordinates loaded - total number of nodes: " << Info::sumValue(mesh.coordinates().clusterSize());

		TimeEvent tElements("elements"); tElements.start();
		elements(mesh._elements);
		materials(mesh._materials);
		tElements.end(); measurement.addEvent(tElements);
		ESINFO(OVERVIEW) << "Elements loaded - total number of elements: " << Info::sumValue(mesh.elements().size());

		mesh.fillNodesFromElements();
		mesh.fillParentElementsToNodes();

		TimeEvent tFaces("faces"); tFaces.start();
		faces(mesh._faces);
		tFaces.end(); measurement.addEvent(tFaces);
		ESINFO(DETAILS) << "Faces loaded - total number of faces: " << Info::sumValue(mesh._faces.size());

		TimeEvent tSettings("settings"); tSettings.start();
		settings(mesh._evaluators, mesh._elements, mesh._faces, mesh._edges, mesh._nodes);
		tSettings.end(); measurement.addEvent(tSettings);

		mesh.fillFacesParents();
		mesh.fillEdgesParents();

		TimeEvent tClusterBoundaries("cluster boundaries"); tClusterBoundaries.start();
		clusterBoundaries(mesh._nodes, mesh._neighbours);
		tClusterBoundaries.end(); measurement.addEvent(tClusterBoundaries);
		ESINFO(OVERVIEW) << "Neighbours loaded - number of neighbours for each cluster is " << Info::averageValue(mesh.neighbours().size());

		close();

		std::vector<std::vector<eslocal> > fPoints;
		fixPoints(fPoints);
		mesh._fixPoints.resize(fPoints.size());
		for (size_t p = 0; p < fPoints.size(); p++) {
			mesh._fixPoints[p].reserve(fPoints[p].size());
			for (size_t i = 0; i < fPoints[p].size(); i++) {
				mesh._fixPoints[p].push_back(mesh.nodes()[fPoints[p][i]]);
			}
		}

		std::vector<eslocal> cornerPoints;
		corners(cornerPoints);
		mesh._corners.reserve(cornerPoints.size());
		for (size_t i = 0; i < cornerPoints.size(); i++) {
			mesh._corners.push_back(mesh.nodes()[cornerPoints[i]]);
		}

		TimeEvent tPartition("partition"); tPartition.start();
		if (partitiate(mesh._partPtrs)) { // manual partition -> map elements to the domains
			mesh.mapElementsToDomains();
			mesh.mapFacesToDomains();
			mesh.mapEdgesToDomains();
			mesh.mapNodesToDomains();
			mesh.mapCoordinatesToDomains();
		}

		if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
			ESTEST(MANDATORY) << "Do not use HYBRID FETI for clusters with 1 domain." << (mesh.parts() > 1 ? TEST_PASSED : TEST_FAILED);
		}

		tPartition.end(); measurement.addEvent(tPartition);
		ESINFO(OVERVIEW) << "Mesh partitioned into " << config::env::MPIsize << " * " << mesh.parts() << " = " << mesh.parts() * config::env::MPIsize
				<< " parts. There is " << intervalStats(mesh._partPtrs) << " elements in subdomain.";

		measurement.totalTime.endWithBarrier(); measurement.printStatsMPI();
	}

	virtual void points(Coordinates &coordinates) = 0;
	virtual void elements(std::vector<Element*> &elements) { }; // Generator, Workbench
	virtual void faces(std::vector<Element*> &faces) { }; // OpenFOAM
	virtual void materials(std::vector<Material> &materials) = 0;
	virtual void clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours) = 0;
	virtual void settings(
			std::vector<Evaluator*> &evaluators,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes) {};

	virtual void open() {};
	virtual void close() {};

	virtual bool partitiate(std::vector<eslocal> &parts)
	{
		mesh.partitiate(config::mesh::SUBDOMAINS);
		return false;
	}

	virtual void fixPoints(std::vector<std::vector<eslocal> > &fixPoints) {};
	virtual void corners(std::vector<eslocal> &corners) {};

protected:
	Loader(Mesh &mesh): mesh(mesh) {};
	virtual ~Loader() {};
	Mesh &mesh;
};

}
}



#endif /* INPUT_LOADER_H_ */
