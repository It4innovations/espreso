
#ifndef INPUT_MESHGENERATOR_FACTORY_H_
#define INPUT_MESHGENERATOR_FACTORY_H_

#include "settings.h"
#include "../generator.h"
#include "../uniformmesh/cube/generator.h"
#include "../uniformmesh/sphere/generator.h"

namespace esinput {

class MeshFactory {

public:
	static Generator* create(int argc, char** argv, size_t index, size_t size);
};

}


#endif /* INPUT_MESHGENERATOR_FACTORY_H_ */
