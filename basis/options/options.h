
#ifndef APP_OPTIONS_OPTIONS_H_
#define APP_OPTIONS_OPTIONS_H_

#include <vector>
#include <iostream>
#include <getopt.h>

struct Options {

	friend std::ostream& operator<<(std::ostream& os, const Options &options);

	Options(int *argc, char*** argv);

	std::string input;
	std::string path;
	size_t verbosity;
	std::vector<std::string> nameless;
};

#endif /* APP_OPTIONS_OPTIONS_H_ */
