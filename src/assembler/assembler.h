
#ifndef SRC_ASSEMBLER_ASSEMBLER_H_
#define SRC_ASSEMBLER_ASSEMBLER_H_

namespace espreso {

struct GlobalConfiguration;
class OldInstance;
class Mesh;

struct Assembler {
	static void compose(const GlobalConfiguration &configuration, OldInstance* &instance, Mesh &mesh);
};

}



#endif /* SRC_ASSEMBLER_ASSEMBLER_H_ */
