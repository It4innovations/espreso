
#include "loader.h"

#include "../mesh/settings/property.h"
#include "../config/description.h"

using namespace espreso::input;

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
	neighbours(mesh._nodes, mesh._neighbours);
	tClusterBoundaries.end(); measurement.addEvent(tClusterBoundaries);
	ESINFO(OVERVIEW) << "Neighbours loaded - number of neighbours for each cluster is " << Info::averageValue(mesh.neighbours().size());


	// LOAD BOUNDARY CONDITIONS
	TimeEvent tSettings("settings"); tSettings.start();
	regions(mesh._evaluators, mesh._regions, mesh._elements, mesh._faces, mesh._edges, mesh._nodes);

	mesh.fillFacesParents();
	mesh.fillEdgesParents();

	boundaryConditions();
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

	if (config::solver::FETI_METHOD == config::solver::FETI_METHODalternative::HYBRID_FETI) {
		ESTEST(MANDATORY) << "Do not use HYBRID FETI for clusters with 1 domain." << (mesh.parts() > 1 ? TEST_PASSED : TEST_FAILED);
	}

	tPartition.end(); measurement.addEvent(tPartition);
	ESINFO(OVERVIEW) << "Mesh partitioned into " << config::env::MPIsize << " * " << mesh.parts() << " = " << mesh.parts() * config::env::MPIsize
			<< " parts. There is " << intervalStats(mesh._partPtrs) << " elements in subdomain.";

	measurement.totalTime.endWithBarrier(); measurement.printStatsMPI();
}

void Loader::boundaryConditions()
{
	auto getValueIndex = [] (const std::vector<std::string> &values, const std::string &parameter) -> size_t {
		if (values.size() == 1 && !Parser::contains(values[0], ":=")) {
			return 0;
		}
		for (size_t i = 0; i < values.size(); i++) {
			if (StringCompare::caseInsensitiveEq(parameter, Parser::strip(Parser::split(values[i], ":=")[0]))) {
				return i;
			}
		}
		return values.size();
	};

	auto loadProperty = [&] (const std::map<std::string, std::string> &regions, const std::vector<std::string> &parameters, const std::vector<Property> &properties) {
		for (auto it = regions.begin(); it != regions.end(); ++it) {
			Region &region = mesh.region(it->first);
			std::vector<std::string> values = Parser::split(it->second, ",;");

			for (size_t p = 0; p < properties.size(); p++) {
				size_t index = getValueIndex(values, parameters[p]);
				if (index < values.size()) {
					std::string value = Parser::contains(values[index], ":=") ? Parser::split(values[index], ":=")[1] : values[index];
					if (value.find("xyzt") == std::string::npos) {
						espreso::Expression expr(value, {});
						mesh._evaluators.push_back(new espreso::ConstEvaluator(expr.evaluate({}), properties[p]));
					} else {
						mesh._evaluators.push_back(new espreso::CoordinatesEvaluator(value, mesh.coordinates(), properties[p]));
					}
					for (size_t i = 0; i < region.elements.size(); i++) {
						region.elements[i]->addSettings(properties[p], mesh._evaluators.back());
					}
				}
			}
		}
	};

	loadProperty(configuration.displacement.values       , { "T", "x", "y", "z" }, { Property::TEMPERATURE, Property::DISPLACEMENT_X, Property::DISPLACEMENT_Y, Property::DISPLACEMENT_Z });
	loadProperty(configuration.normal_presure.values     , { "P" }               , { Property::PRESSURE });
	loadProperty(configuration.translation_motions.values, { "x", "y", "z" }     , { Property::TRANSLATION_MOTION_X, Property::TRANSLATION_MOTION_Y, Property::TRANSLATION_MOTION_Z });
	loadProperty(configuration.acceleration.values       , { "x", "y", "z" }     , { Property::ACCELERATION_X, Property::ACCELERATION_Y, Property::ACCELERATION_Z });
	loadProperty(configuration.initial_temperature.values, { }                   , { Property::INITIAL_TEMPERATURE });
	loadProperty(configuration.temperature.values        , { }                   , { Property::TEMPERATURE });
	loadProperty(configuration.heat_source.values        , { "T" }               , { Property::HEAT_SOURCE });
	loadProperty(configuration.thickness.values          , { }                   , { Property::THICKNESS });
	loadProperty(configuration.obstacle.values           , { }                   , { Property::OBSTACLE });
	loadProperty(configuration.normal_direction.values   , { }                   , { Property::NORMAL_DIRECTION });

	for (auto it = configuration.material_set.values.begin(); it != configuration.material_set.values.end(); ++it) {
		Region &region = mesh.region(it->second);
		for (size_t e = 0; e << region.elements.size(); e++) {
			region.elements[e]->setParam(Element::MATERIAL, it->first -1);
		}
		for (auto p = configuration.materials.configurations[it->first - 1]->parameters.begin(); p != configuration.materials.configurations[it->first - 1]->parameters.end(); ++p) {
			mesh._materials[it->first - 1].setParameter(p->first, p->second->get());
		}
	}

}



