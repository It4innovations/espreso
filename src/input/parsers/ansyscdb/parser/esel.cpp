
#include "esel.h"

#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

ESel::ESel()
: type(Type::UNKNOWN), item(Item::UNKNOWN), comp(Comp::UNKNOWN),
  VMIN(0), VMAX(0), VINC(1),
  KABS(false)
{

}

ESel& ESel::parse(const char* begin)
{
	std::string commandLine = Parser::removecomments(Parser::getLine(begin), "!");
	std::vector<std::string> command = Parser::split(Parser::strip(commandLine), ",", false);

	switch (command.size()) {
	case 8:
		KABS = std::stoi(command[7]);
	case 7:
		if (command[6].size()) {
			VINC = std::stoi(command[6]);
		}
	case 6:
		if (command[5].size()) {
			size_t end;
			VMAX = std::stoi(command[5], &end) - 1;
			if (end != command[5].size()) {
				eslog::error("ESPRESO Workbench parser error: implement esel, VMAX='%s'\n", command[5].c_str());
			}
		}
	case 5:
		if (command[4].size()) {
			size_t end;
			VMIN = std::stoi(command[4], &end) - 1;
			if (end != command[4].size()) {
				eslog::error("ESPRESO Workbench parser error: implement esel, VMIN='%s'\n", command[4].c_str());
			}
		}
		if (command.size() < 6 || command[5].size() == 0) {
			VMAX = VMIN;
		}
	case 4:
		if (command[3].size()) {
			eslog::error("ESPRESO Workbench parser error: implement esel, Comp='%s'\n", command[3].c_str());
		}
	case 3:
		if (StringCompare::caseInsensitiveEq("ELEM", command[2])) {
			item = Item::ELEM;
		}
		if (StringCompare::caseInsensitiveEq("ADJ", command[2])) {
			item = Item::ADJ;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("CENT", command[2])) {
			item = Item::CENT;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("TYPE", command[2])) {
			item = Item::TYPE;
		}
		if (StringCompare::caseInsensitiveEq("ENAME", command[2])) {
			item = Item::ENAME;
			eslog::error("ESPRESO Workbench parser error: not implemented esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("MAT", command[2])) {
			item = Item::MAT;
		}
		if (StringCompare::caseInsensitiveEq("REAL", command[2])) {
			item = Item::REAL;
		}
		if (StringCompare::caseInsensitiveEq("ESYS", command[2])) {
			item = Item::ESYS;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("PART", command[2])) {
			item = Item::PART;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("LIVE", command[2])) {
			item = Item::LIVE;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("LAYER", command[2])) {
			item = Item::LAYER;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("SEC", command[2])) {
			item = Item::SEC;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("STRA", command[2])) {
			item = Item::STRA;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("SFE", command[2])) {
			item = Item::SFE;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("BFE", command[2])) {
			item = Item::BFE;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("PATH", command[2])) {
			item = Item::PATH;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
		if (StringCompare::caseInsensitiveEq("ETAB", command[2])) {
			item = Item::ETAB;
			eslog::error("ESPRESO Workbench parser error: implement esel, Item='%s'\n", command[2].c_str());
		}
	case 2:
		if (StringCompare::caseInsensitiveEq("S", command[1])) {
			type = Type::S;
		}
		if (StringCompare::caseInsensitiveEq("R", command[1])) {
			type = Type::R;
			eslog::error("ESPRESO Workbench parser error: implement esel, Type='%s'\n", command[1].c_str());
		}
		if (StringCompare::caseInsensitiveEq("A", command[1])) {
			type = Type::A;
		}
		if (StringCompare::caseInsensitiveEq("U", command[1])) {
			type = Type::U;
			eslog::error("ESPRESO Workbench parser error: implement esel, Type='%s'\n", command[1].c_str());
		}
		if (StringCompare::caseInsensitiveEq("ALL", command[1])) {
			type = Type::ALL;
		}
		if (StringCompare::caseInsensitiveEq("NONE", command[1])) {
			type = Type::NONE;
		}
		if (StringCompare::caseInsensitiveEq("INVE", command[1])) {
			type = Type::INVE;
			eslog::error("ESPRESO Workbench parser error: implement esel, Type='%s'\n", command[1].c_str());
		}
		if (StringCompare::caseInsensitiveEq("STAT", command[1])) {
			type = Type::STAT;
			eslog::error("ESPRESO Workbench parser error: implement esel, Type='%s'\n", command[1].c_str());
		}
		break;
	default:
		eslog::error("ESPRESO Workbench parser error: unknown format of '%s'\n", commandLine.c_str());
	}

	if (type == Type::UNKNOWN) {
		eslog::error("ESPRESO Workbench parser error: implement esel, Type='%s'\n", command[1].c_str());
	}
	if (comp != Comp::UNKNOWN) {
		eslog::error("ESPRESO Workbench parser error: implement esel, Comp='%s'\n", command[3].c_str());
	}

	WorkbenchParser::fillIndices(begin, begin, begin);
	return *this;
}


