
#include "settings.h"

using namespace espreso::input;

void CubesSettings::defaultCubesSettings()
{


}

CubesSettings::CubesSettings(const ArgsConfiguration &configuration, size_t index, size_t size)
: cube{CubeSettings(index, size, "MESH1_"), CubeSettings(index, size, "MESH2_")}
{
	parameters.insert(parameters.end(), cube[0].parameters.begin(), cube[0].parameters.end());
	parameters.insert(parameters.end(), cube[1].parameters.begin(), cube[1].parameters.end());
	defaultCubesSettings();
	ESINFO(OVERVIEW) << "Load cubes setting from file " << configuration.path;
	ParametersReader::fromConfigurationFileWOcheck(configuration, parameters);
}

CubesSettings::CubesSettings(size_t index, size_t size)
: cube{CubeSettings(index, size, "MESH1_"), CubeSettings(index, size, "MESH2_")}
{
	cube[0] = CubeSettings(index, size);
	cube[1] = CubeSettings(index, size);
	parameters.insert(parameters.end(), cube[0].parameters.begin(), cube[0].parameters.end());
	parameters.insert(parameters.end(), cube[1].parameters.begin(), cube[1].parameters.end());
	defaultCubesSettings();
}




