
#ifndef INPUT_ANSYS_WORKBENCH_PARSER_H_
#define INPUT_ANSYS_WORKBENCH_PARSER_H_

#include <map>

#include "esbasis.h"
#include "esmesh.h"
#include "../utils.h"

namespace espreso {
namespace input {

enum class WorkbenchCommands {
	END,
	WB,
	NBLOCK,
	EBLOCK,
	MP,
	DISPLACEMENT,
	FORCE,
	LOADVAR,
	DIM,

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

	void open(std::string path)
	{
		_file.open(path.c_str());
		if (!_file.is_open()) {
			ESINFO(GLOBAL_ERROR) << "Cannot load mesh from file: " << path;
		}
	}

	void close()
	{
		_file.close();
	}

	WorkbenchCommands process(); // process to a next command that returns some data

	bool workbench(const std::string type, const std::string status);

	void nblock(Coordinates &coordinates);

	void eblock(std::vector<Element*> &elements);
	void mp(std::vector<Material> &materials);

	void eblock(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void cmblock(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void displacement(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void force(std::vector<Evaluator*> &evaluators, std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges, std::vector<Element*> &nodes);
	void loadvar();

	~WorkbenchParser()
	{
		for (size_t e = 0; e < _evaluators.size(); e++) {
			delete _evaluators[e];
		}
		for (size_t t = 0; t < _tables.size(); t++) {
			delete _tables[t];
		}
	}

protected:
	void et();
	void cmsel();
	void nsel();
	void esel();
	void dim();

	std::vector<std::string> divide(std::string &line, std::string delim = ",");
	std::vector<int> parseBlockHeader(std::string &line);
	bool trash(const std::string &line);

	std::ifstream _file;
	std::string _line;

	std::map<std::string, WorkbenchCommands, CaseInsensitiveCompare> _commands;

	int bodyCounter;
	std::vector<int> eType;
	std::vector<std::pair<std::string, ConditionElements> > selections;
	Mesh &_mesh;
	std::vector<Evaluator*> _evaluators;
	std::map<std::string, std::vector<Element*> > _regions;
	std::vector<TableEvaluator*> _tables;
};

}
}

#endif /* INPUT_ANSYS_WORKBENCH_PARSER_H_ */
