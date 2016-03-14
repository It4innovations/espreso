
#ifndef APP_OPTIONS_OPTIONS_H_
#define APP_OPTIONS_OPTIONS_H_

#include "mpi.h"

#include <vector>
#include <iostream>
#include <getopt.h>

#include "esconfig.h"
#include "../logging/logging.h"


namespace espreso {

struct Options {

	friend std::ostream& operator<<(std::ostream& os, const Options &options);

	Options(int *argc, char*** argv);

	void configure();

	std::string input;
	std::string path;
	size_t verboseLevel;
	size_t testingLevel;
	std::vector<std::string> nameless;
};

}

#endif /* APP_OPTIONS_OPTIONS_H_ */
