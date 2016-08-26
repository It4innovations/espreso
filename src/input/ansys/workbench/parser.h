
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
	LOADVAR,

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

//	void eblock(std::vector<Dirichlet*> &dirichlet);
//	void cmblock(std::vector<Dirichlet*> &dirichlet);
//	void displacement(std::vector<Dirichlet*> &dirichlet);
	void loadvar();

protected:
	void et();
	void cmsel();
	void nsel();
	void esel();

	std::vector<std::string> divide(std::string &line, std::string delim = ",");
	std::vector<int> parseBlockHeader(std::string &line);
	bool trash(const std::string &line);

	std::ifstream _file;
	std::string _line;

	std::map<std::string, WorkbenchCommands, CaseInsensitiveCompare> _commands;

	int bodyCounter;
	std::vector<int> eType;
	std::vector<std::pair<std::string, ConditionElements> > selections;
	int nSelection;
	int eSelection;
	Mesh &_mesh;
};

}
}

#endif /* INPUT_ANSYS_WORKBENCH_PARSER_H_ */
