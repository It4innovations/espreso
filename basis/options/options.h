
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

	Options(): verboseLevel(0), testingLevel(0), measureLevel(0) {};
	Options(int *argc, char*** argv);

	void configure();
	void setFromFile(const std::string &file);

	std::string input;
	std::string path;
	size_t verboseLevel;
	size_t testingLevel;
	size_t measureLevel;
	std::vector<std::string> nameless;
};

}

#endif /* APP_OPTIONS_OPTIONS_H_ */
