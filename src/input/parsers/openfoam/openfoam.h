
#ifndef SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_
#define SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_

#include "input/input.h"
#include "parser/foamfileheader.h"

namespace espreso {

struct NamedData;
class InputConfiguration;

class InputOpenFoam: public Input {
public:
	InputOpenFoam();
	~InputOpenFoam();

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);

	void initVariables(Mesh &mesh);
	void finishVariables();
	int timeSteps();
	void nextTimeStep(Mesh &mesh);

protected:
	std::vector<int> nvariables, svariables, evariables;
	std::vector<FoamFileHeader> vheader;
	AsyncFilePack variablePack;

	std::string variablePath;
	size_t timestep;
	std::vector<std::string> timesteps, variableNames;

	InputOpenFoam *loader;
};

class InputOpenFoamSequential: public InputOpenFoam {
public:
	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);

	void initVariables(Mesh &mesh);
	void finishVariables();
	void nextTimeStep(Mesh &mesh);
	void ivariables(const InputConfiguration &configuration);

protected:
	NodesBlocks nodes;
	FacesBlocks faces;
	OrderedRegions regions;
};

class InputOpenFoamParallel: public InputOpenFoam {
public:
	InputOpenFoamParallel(esint domains): domains(domains) {}
	esint domains;

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);

	void initVariables(Mesh &mesh);
	void finishVariables();
	void nextTimeStep(Mesh &mesh);
	void ivariables(const InputConfiguration &configuration);

protected:
	NodesBlocks nodes;
	FacesBlocks faces;
	VariablesBlocks variables;
	OrderedRegions regions;
};

class InputOpenFoamParallelDirect: public InputOpenFoam {
public:
	InputOpenFoamParallelDirect(esint domains): domains(domains) {}
	esint domains;

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);

	void initVariables(Mesh &mesh);
	void finishVariables();
	void nextTimeStep(Mesh &mesh);
	void ivariables(const InputConfiguration &configuration);

protected:
	NodesDomain nodes;
	Faces faces;
	OrderedRegions regions;
};

}

#endif /* SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_ */
