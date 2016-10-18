
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

void AnsysWorkbench::settings(
			std::vector<Evaluator*> &evaluators,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes)
{
	while (true) {
		switch (_parser.process()) {
		case WorkbenchCommands::WB:
			if (_parser.workbench("load", "end")) {
				// skip end because there can be another settings
				// return;
			}
			break;
		case WorkbenchCommands::CMBLOCK:
			_parser.cmblock(elements, faces, edges, nodes);
			break;
		case WorkbenchCommands::DISPLACEMENT:
			_parser.displacement(evaluators, elements, faces, edges, nodes);
			break;
		case WorkbenchCommands::FORCE:
			_parser.force(evaluators, elements, faces, edges, nodes);
			break;
		case WorkbenchCommands::OBSTACLE:
			_parser.obstacle(evaluators, elements, faces, edges, nodes);
			break;
		case WorkbenchCommands::EBLOCK:
			_parser.eblock(elements, faces, edges, nodes);
			break;
		default:
			return;
		}
	}
}


void AnsysWorkbench::clusterBoundaries(std::vector<Element*> &nodes, std::vector<int> &neighbours)
{
	for (size_t i = 0; i < mesh.coordinates().clusterSize(); i++) {
		nodes[i]->clusters().push_back(0);
	}
}


