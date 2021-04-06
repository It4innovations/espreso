
#ifndef SRC_INPUT_ABAQUS_ABAQUS_H_
#define SRC_INPUT_ABAQUS_ABAQUS_H_

#include "input/meshbuilder.h"
#include "basis/io/inputfile.h"

#include <cstddef>
#include <vector>

namespace espreso {

class InputConfiguration;
class Mesh;

struct NList;
struct EList;
struct CMBlock;
struct ET;
struct Eset;
struct Nset;
struct CM;
struct BlockFinish;
struct SSection;
struct AbaqusMaterial;
struct Elemat;

class AbaqusLoader: public MeshBuilder {
public:
	AbaqusLoader(const InputConfiguration &configuration);
	void load();

protected:
	void scan();
	void parse();

	const InputConfiguration &_configuration;

	std::vector<NList> _NLists;
	std::vector<EList> _ELists;
	std::vector<BlockFinish> _blockFinishs;
	std::vector<Eset> _Esets;
	std::vector<SSection> _SSections;
	std::vector<AbaqusMaterial> _Materials;
	std::vector<Elemat> _Elemats;
	std::vector<Nset> _Nsets;
	InputFilePack _pfile;
};

}



#endif /* SRC_INPUT_ABAQUS_ABAQUS_H_ */
