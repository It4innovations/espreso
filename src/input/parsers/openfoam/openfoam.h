
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
	double nextVariables(Mesh &mesh);

	int variables();

protected:
	std::vector<DatabaseOffset> nblocks, eblocks;

	std::vector<esint> ndistribution, edistribution;
	std::vector<int> nvariables, svariables, evariables;
	std::vector<FoamFileHeader> vheader;
	AsyncFilePack variablePack;

	std::string variablePath;
	size_t timestep;
	std::vector<std::string> timesteps, variableNames;
	std::vector<esint> variableRequests, variableRequestsMap;
	std::vector<NamedData*> vdata;

	InputOpenFoam *loader;
};

class InputOpenFoamSequential: public InputOpenFoam {
public:
	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);
	double nextVariables(Mesh &mesh);
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
	double nextVariables(Mesh &mesh);
	void ivariables(const InputConfiguration &configuration);

protected:
	NodesBlocks nodes;
	FacesBlocks faces;
	OrderedRegions regions;

	ivector<esint> nperm, eperm;
};

class InputOpenFoamParallelDirect: public InputOpenFoam {
public:
	InputOpenFoamParallelDirect(esint domains): domains(domains) {}
	esint domains;

	void load(const InputConfiguration &configuration);
	void build(Mesh &mesh);
	double nextVariables(Mesh &mesh);
	void ivariables(const InputConfiguration &configuration);

protected:
	NodesDomain nodes;
	Faces faces;
	OrderedRegions regions;
};

}

#endif /* SRC_INPUT_FORMATS_OPENFOAM_OPENFOAM_H_ */
