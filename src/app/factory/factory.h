
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include "esmesh.h"
#include "esassembler.h"

namespace espreso {

struct GlobalConfiguration;

struct Factory {

	Factory(const GlobalConfiguration &configuration);
	~Factory()
	{
		delete instance;
	}

	void solve(const std::string &outputFile);

	double norm() const;

	Instance *instance;
	Mesh mesh;

private:
	std::vector<std::vector<double> > _solution;
};

}



#endif /* APP_FACTORY_FACTORY_H_ */
