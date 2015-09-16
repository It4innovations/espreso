
#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include "parameter.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>

class Configuration {

public:

	enum InputParameter {
		// General settings
		CMD_LINE_ARGUMENTS,
		// Ansys settings
		ANSYS_DIR,
		ANSYS_DIRICHLET_X,
		ANSYS_DIRICHLET_Y,
		ANSYS_DIRICHLET_Z,
		ANSYS_FORCES_X,
		ANSYS_FORCES_Y,
		ANSYS_FORCES_Z,
		// PermonCube settings
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
		// Solver settings
		SOLVER_ITERATIONS,

		// Dummy attribute -> attributes counter
		ATTRIBUTES_COUNT
	};

	Configuration(std::string configFile);
	Configuration(): _parameters(ATTRIBUTES_COUNT)
	{
		_fillDefaultValues();
	}

	void print();
	void description();

	template<class ParameterType>
	ParameterType value(InputParameter parameter);

	template<class ParameterType>
	ParameterType value(size_t parameter)
	{
		return value<ParameterType>(static_cast<InputParameter>(parameter));
	}

private:

	void _fillDefaultValues();

	std::vector<Parameter*> _parameters;

};


#endif /* CONFIGURATION_H_ */
