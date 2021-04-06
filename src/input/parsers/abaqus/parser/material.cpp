
#include "material.h"

#include "basis/containers/point.h"
#include "basis/containers/tarray.h"
#include "basis/utilities/parser.h"

using namespace espreso;

size_t AbaqusMaterial::size = 9;
const char* AbaqusMaterial::upper     = "*MATERIAL";
const char* AbaqusMaterial::lower     = "*material";
const char* AbaqusMaterial::sentence  = "*Material";
AbaqusMaterial::AbaqusMaterial()
:youngs_modulus(1),poisson_ratio(0.3)
{
	memset(Name, '\0', MAX_NAME_SIZE);

}

AbaqusMaterial& AbaqusMaterial::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);

	std::vector<std::string> command = Parser::split(commandLine, ",", false);
	for (size_t i = 0; i < command.size(); i++)
	{
		std::vector<std::string> material_line = Parser::split(command[i], "=", false);
		material_line[0] = Parser::strip(material_line[0]);

		if ( StringCompare::caseInsensitiveEq(material_line[0], "NAME")) {
			material_line[1] = Parser::strip(material_line[1]);
			memcpy(Name,material_line[1].data(), material_line[1].size());
		}

		const char* linebegin = begin ;
		if (*(linebegin ) != '\n') {
			while (*linebegin++ != '\n'); // start at new line
		}
		if (*(linebegin ) != '\n') {
			while (*linebegin++ != '\n'); // start at new line
		}
		std::string propertyLine = Parser::getLine(linebegin);
		std::vector<std::string> property = Parser::split(propertyLine, ",", false);
		youngs_modulus = atof(property[0].data());
		poisson_ratio = atof(property[1].data());
	}

	return *this;
}
