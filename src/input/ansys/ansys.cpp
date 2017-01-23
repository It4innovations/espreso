
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


void AnsysWorkbench::elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges)
{
	while (true) {
		switch (_parser.process()) {
		case WorkbenchCommands::WB:
			if (_parser.workbench("elem", "end")) {
				return;
			}
			break;
		case WorkbenchCommands::EBLOCK: {
			std::vector<Region*> dummyRegions;
			_parser.eblock(elements, dummyRegions, faces, edges);
			break;
		}
		case WorkbenchCommands::END:
			return;
		default:
			break;
		}
	}
}

void AnsysWorkbench::materials(std::vector<Material*> &materials)
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
		default:
			break;
		}
	}
}

void AnsysWorkbench::regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region*> &regions,
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
			_parser.cmblock(elements, regions, faces, edges, nodes);
			break;
		case WorkbenchCommands::DIRICHLET:
			_parser.dirichlet(evaluators, regions, elements, faces, edges, nodes);
			break;
		case WorkbenchCommands::FORCE:
			_parser.force(evaluators, regions, elements, faces, edges, nodes);
			break;
		case WorkbenchCommands::SURFACE_EFFECT:
			_parser.sf(evaluators, regions, elements, faces, edges);
			break;
		case WorkbenchCommands::ACCELERATION:
			_parser.acceleration(evaluators, regions);
			break;
		case WorkbenchCommands::INITIAL_TEMPERATURE:
			_parser.initial_temperature(evaluators, regions);
			break;
		case WorkbenchCommands::OBSTACLE:
			_parser.obstacle(evaluators, regions, elements, faces, edges, nodes);
			break;
		case WorkbenchCommands::EBLOCK:
			_parser.eblock(elements, regions, faces, edges);
			break;
		default:
			return;
		}
	}
}


void AnsysWorkbench::neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours, const std::vector<Element*> &faces, const std::vector<Element*> &edges)
{
	for (size_t i = 0; i < mesh.coordinates().clusterSize(); i++) {
		nodes[i]->clusters().push_back(0);
	}
}


