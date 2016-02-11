
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include "esmesh.h"
#include "esinput.h"
#include "esoutput.h"
#include "esassemblers.h"

class Factory {

public:
	Factory(const Options &options);

	void solve( eslocal steps );
	void store(const std::string &file);

	~Factory();

private:
	assembler::AssemblerBase *_assembler;
	std::vector<std::vector<double> > _solution;

	mesh::Mesh *_mesh;
	mesh::Mesh *_surface;
	assembler::APIHolder *_apiHolder;
};



#endif /* APP_FACTORY_FACTORY_H_ */
