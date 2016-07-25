
#include "ansys.h"
using namespace espreso::input;

void AnsysWorkbench::points(Coordinates &coordinates)
{
	while (true) {
		switch (_parser.process()) {
		case WorkbenchCommands::WB:
			break;
		case WorkbenchCommands::NBLOCK: {
			_parser.nblock(coordinates);
			return;
		}
		default:
			return;
		}
	}
}


void AnsysWorkbench::elements(std::vector<Element*> &elements)
{
	while (true) {
		switch (_parser.process()) {
		case WorkbenchCommands::WB:
			if (_parser.workbench("elem", "end")) {
				return;
			}
			break;
		case WorkbenchCommands::EBLOCK: {
			_parser.eblock(elements);
			break;
		}
		case WorkbenchCommands::END:
			return;
		}
	}
}

void AnsysWorkbench::materials(std::vector<Material> &materials)
{
	while (true) {
		switch (_parser.process()) {
		case WorkbenchCommands::WB:
			if (_parser.workbench("mat", "end")) {
				return;
			}
			break;
		case WorkbenchCommands::MP: {
			_parser.mp(materials);
			break;
		}
		case WorkbenchCommands::END:
			return;
		}
	}
}

void AnsysWorkbench::boundaryConditions(Coordinates &coordinates, std::vector<BoundaryCondition*> &conditions)
{
	while (true) {
		switch (_parser.process()) {
		case WorkbenchCommands::WB:
			if (_parser.workbench("load", "end")) {
				return;
			}
			break;
		case WorkbenchCommands::CMBLOCK:
			_parser.cmblock(conditions);
			break;
		case WorkbenchCommands::DISPLACEMENT:
			_parser.displacement(conditions);
			break;
		case WorkbenchCommands::EBLOCK:
			_parser.eblock(conditions);
			break;
		default:
			return;
		}
	}
}


void AnsysWorkbench::clusterBoundaries(Boundaries &boundaries, std::vector<int> &neighbours)
{
	boundaries.resize(mesh.coordinates().clusterSize());
	for (size_t i = 0; i < mesh.coordinates().clusterSize(); i++) {
		boundaries[i].push_back(0);
	}
}


