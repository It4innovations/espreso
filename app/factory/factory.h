
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include "esmesh.h"
#include "esinput.h"
#include "esoutput.h"
#include "esassemblers.h"

class Factory {

public:
	Factory(int *argc, char ***argv);

	void solve( eslocal steps );

	~Factory();

private:
	assembler::AssemblerBase *_assembler;
	mesh::Mesh *_mesh;
	mesh::SurfaceMesh *_surface;
};



#endif /* APP_FACTORY_FACTORY_H_ */
