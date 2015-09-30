
#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "parameter.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>

class Configuration {

public:

	enum ElementType {
		HEXA8,
		HEXA20,
		TETRA4,
		TETRA10,
		PRISMA6,
		PRISMA15,
		PYRAMID5,
		PYRAMID13
	};

	enum PermonCubeShapes {
		CUBE,
		SPHERE
	};

	enum InputParameter {
		// General settings
		CMD_LINE_ARGUMENTS,

		// Mesh settings
		MESH_FILE,
		BOUNDARIES_FILE,
		MESH_SUBDOMAINS,
		MESH_FIX_POINTS,
		MESH_CORNERS_NUMBER,
		MESH_CORNERS_IN_CORNER,
		MESH_CORNERS_IN_EDGES,
		MESH_CORNERS_IN_FACES,

		// Ansys settings
		ANSYS_DIR,
		ANSYS_DIRICHLET_X,
		ANSYS_DIRICHLET_Y,
		ANSYS_DIRICHLET_Z,
		ANSYS_FORCES_X,
		ANSYS_FORCES_Y,
		ANSYS_FORCES_Z,

		// PermonCube settings
		PMCUBE_SHAPE,
		PMCUBE_ELEMENT_TYPE,
		PMCUBE_CLUSTERS_X,
		PMCUBE_CLUSTERS_Y,
		PMCUBE_CLUSTERS_Z,
		PMCUBE_SUBDOMAINS_X,
		PMCUBE_SUBDOMAINS_Y,
		PMCUBE_SUBDOMAINS_Z,
		PMCUBE_ELEMENTS_X,
		PMCUBE_ELEMENTS_Y,
		PMCUBE_ELEMENTS_Z,
		PMCUBE_FIX_ZERO_PLANES,
		PMCUBE_FIX_BOTTOM,
		PMCUBE_CORNERS_X,
		PMCUBE_CORNERS_Y,
		PMCUBE_CORNERS_Z,
		PMCUBE_SPHERE_LAYERS,
		PMCUBE_SPHERE_INNER_RADIUS,
		PMCUBE_SPHERE_OUTER_RADIUS,

		// Solver settings
		SOLVER_ITERATIONS,

		// Dummy attribute -> attributes counter
		ATTRIBUTES_COUNT
	};

	Configuration(std::string configFile, int argc, char** argv);
	Configuration(): _parameters(ATTRIBUTES_COUNT)
	{
		fillDefaultValues();
		configure(0, NULL);
	}
	Configuration(const Configuration &other)
	{
		_parameters.resize(other._parameters.size());
		for (size_t i = 0; i < other._parameters.size(); i++) {
			_parameters[i] = other._parameters[i]->copy();
		}
	}
	~Configuration();

	void print() const;
	void description() const;

	Parameter* parameter(InputParameter parameter)
	{
		return _parameters[parameter];
	}

	Parameter* parameter(size_t parameter)
	{
		return _parameters[parameter];
	}

	template<class ParameterType>
	const ParameterType& value(InputParameter parameter) const;

	template<class ParameterType>
	const ParameterType& value(size_t parameter) const
	{
		return value<ParameterType>(static_cast<InputParameter>(parameter));
	}

private:

	void fillDefaultValues();
	void configure(int argc, char** argv);

	std::vector<Parameter*> _parameters;
};


#endif /* CONFIGURATION_H_ */
