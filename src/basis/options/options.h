
#ifndef APP_OPTIONS_OPTIONS_H_
#define APP_OPTIONS_OPTIONS_H_

#include "mpi.h"

#include <vector>
#include <iostream>
#include <getopt.h>
#include <signal.h>
#include <csignal>

#include "esconfig.h"
#include "../logging/logging.h"


namespace espreso {

struct Options {

	friend std::ostream& operator<<(std::ostream& os, const Options &options);

	Options(): verbose(0), measure(0), testing(0), help(0) { };
	Options(int *argc, char*** argv);

	void configure();

	std::string executable;
	std::string path;
	std::string input;
	std::vector<std::string> nameless;

private:
	size_t verbose;
	size_t measure;
	size_t testing;
	size_t help;
};

}

#endif /* APP_OPTIONS_OPTIONS_H_ */
