
#ifndef INPUT_ANSYS_PARSER_H_
#define INPUT_ANSYS_PARSER_H_

#include <map>
#include <vector>
#include <fstream>

#include "utils.h"
#include "../../basis/utilities/parser.h"

namespace espreso {

enum class Property;
class Mesh;
class Coordinates;
class Region;
class Material;
class Evaluator;
class TableEvaluator;

namespace input {

enum class WorkbenchCommands {
	END,
	WB,
	NBLOCK,
	EBLOCK,
	MP,
	DIRICHLET,
	INITIAL_TEMPERATURE,
	FORCE,
	SURFACE_EFFECT,
	ACCELERATION,
	OBSTACLE,
	LOADVAR,
	DIM,
	MPTEMP,

	CMBLOCK,
	ET,
	CMSEL,
	NSEL,
	ESEL
};

enum ConditionElements {
	NODES,
	FACES,
	ELEMENTS
};

class WorkbenchParser {

public:
	WorkbenchParser(Mesh &mesh);

	void open(std::string path);
	void close()
	{
		_file.close();
	}

	WorkbenchCommands process(); // process to a next command that returns some data

	bool workbench(const std::string type, const std::string status);
	void nblock(Coordinates &coordinates);
	void eblock(std::vector<Element*> &elements, std::vector<Region*> &regions, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void mp(std::vector<Material*> &materials, Evaluator *evaluator = NULL);
	void mptemp(std::vector<Material*> &materials);
	void cmblock(std::vector<Element*> &elements, std::vector<Region*> &regions, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void dirichlet(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void force(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void acceleration(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions);
	void initial_temperature(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions);
	void obstacle(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void sf(std::vector<Evaluator*> &evaluators, std::vector<Region*> &regions, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	void loadvar();

	~WorkbenchParser();

protected:
	void et();
	void cmsel();
	void nsel();
	void esel();
	void dim();

	std::vector<std::string> divide(std::string &line, std::string delim = ",");
	std::vector<int> parseBlockHeader(std::string &line);
	bool trash(const std::string &line);

	bool setProperty(const std::string &parameter, const std::string &value, const std::string &name, espreso::Property property, size_t loadStep, espreso::Region *region, std::vector<Evaluator*> &evaluators);

	std::ifstream _file;
	std::string _line;

	std::map<std::string, WorkbenchCommands, StringCompare> _commands;

	int bodyCounter;
	std::vector<int> eType;
	Mesh &_mesh;

	std::string _selectedRegion;
	std::vector<TableEvaluator*> _tables;
};

}
}

#endif /* INPUT_ANSYS_PARSER_H_ */
