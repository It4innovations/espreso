
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include "esmesh.h"
#include "esinput.h"
#include "esoutput.h"
#include "esassemblers.h"

namespace espreso {

class Factory {

public:
	Factory(const Options &options);

	void solve();
	void store(const std::string &file);

	~Factory();

private:
	AssemblerBase *_assembler;
	std::vector<std::vector<double> > _solution;

	Mesh *_mesh;
	Mesh *_surface;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
