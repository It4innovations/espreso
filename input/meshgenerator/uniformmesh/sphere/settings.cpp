
#include "settings.h"

using namespace espreso::input;

static std::vector<Description> createSphereSetting()
{
	std::vector<Description> description(UniformSettings::description);

	description.push_back({
		INTEGER_PARAMETER, "LAYERS", "Number of layers of the sphere."
	});
	description.push_back({
		INTEGER_PARAMETER, "GRID", "Grid size of one side of the sphere."
	});
	description.push_back({
		DOUBLE_PARAMETER, "INNER_RADIUS", "Inner radius of the sphere."
	});
	description.push_back({
		DOUBLE_PARAMETER, "OUTER_RADIUS", "Outer radius of the sphere."
	});

	std::vector<std::pair<std::string, std::string> > axis = {
			{ "X", "x" },
			{ "Y", "y" },
			{ "Z", "z" }
	};
	std::vector<std::pair<std::string, std::string> > properties = {
			{ "DIRICHLET", "Dirichlet" },
			{ "FORCES", "Force" }
	};
	std::vector<std::pair<std::string, std::string> > sphere_faces = {
			{ "INNER", "inner" },
			{ "OUTER", "outer" }
	};

	for (size_t i = 0; i < axis.size(); i++) {
		for (size_t j = 0; j < properties.size(); j++) {
			for (size_t k = 0; k < sphere_faces.size(); k++) {
				description.push_back({
					DOUBLE_PARAMETER,
					properties[j].first + "_" + sphere_faces[k].first + "_" + axis[i].first,
					properties[j].second + " on the " + sphere_faces[k].second + " face in " + axis[i].second + "-axis."
				});
			}
		}
	}

	return description;
};

std::vector<Description> SphereSettings::description = createSphereSetting();

SphereSettings::SphereSettings(const Options &options, size_t index, size_t size)
: UniformSettings(options, index, size)
{
	ESINFO(OVERVIEW) << "Load sphere setting from file " << options.path;
	Configuration configuration(SphereSettings::description, options);

	layers = configuration.value("LAYERS", 1);
	grid = configuration.value("GRID", 1);
	innerRadius = configuration.value("INNER_RADIUS", 9.0);
	outerRadius = configuration.value("OUTER_RADIUS", 12.0);

	std::vector<std::string> axis = { "X", "Y", "Z" };
	std::vector<std::string> properties = { "DIRICHLET", "FORCES" };
	std::vector<std::string> sphere_faces = { "INNER", "OUTER" };

	fillCondition.resize(sphere_faces.size());
	boundaryCondition.resize(sphere_faces.size());

	for (size_t f = 0; f < sphere_faces.size(); f++) {
		for (size_t p = DIRICHLET_X; p <= FORCES_Z; p++) {
			std::string name = properties[p / 3] + "_" + sphere_faces[f] + "_" + axis[p % 3];
			fillCondition[f][p] = configuration.isSet(name);
			boundaryCondition[f][p] = configuration.value<double>(name, 0);
		}
	}
}

SphereSettings::SphereSettings(size_t index, size_t size)
: UniformSettings(index, size)
{
	layers = 1;
	grid = 1;
	innerRadius = 9;
	outerRadius = 12;

	std::vector<std::string> axis = { "X", "Y", "Z" };
	std::vector<std::string> properties = { "DIRICHLET", "FORCES" };
	std::vector<std::string> sphere_faces = { "INNER", "OUTER" };

	fillCondition.resize(sphere_faces.size());
	boundaryCondition.resize(sphere_faces.size());

	for (size_t f = 0; f < sphere_faces.size(); f++) {
		for (size_t p = DIRICHLET_X; p <= FORCES_Z; p++) {
			fillCondition[f][p] = false;
			boundaryCondition[f][p] = 0;
		}
	}
}



