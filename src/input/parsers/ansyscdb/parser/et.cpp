
#include "et.h"

#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

ET::ET()
: id(-1), type(-1)
{

}

ET& ET::parse(const char* begin)
{
	std::string commandLine = Parser::getLine(begin);
	std::vector<std::string> command = Parser::split(Parser::strip(commandLine), ",", false);

	if (command.size() < 3) {
		eslog::error("ESPRESO Workbench parser error: unknown et format of '%s'\n", commandLine.c_str());
	}

	if (
			!StringCompare::caseInsensitiveEq(command[1], "tid") &&
			!StringCompare::caseInsensitiveEq(command[1], "_tid") &&
			!StringCompare::caseInsensitiveEq(command[1], "cid")) {
		id = std::stoi(command[1]) - 1;
	}
	type = std::stoi(command[2]);

	WorkbenchParser::fillIndices(begin, begin, begin);

	if (id != -1 && etype() == ETYPE::UNKNOWN) {
		eslog::error("ESPRESO Workbench parser error: unknown et type '%s'\n", command[2].c_str());
	}

	return *this;
}


