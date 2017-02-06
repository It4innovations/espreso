
#ifndef SRC_CONFIGURATION_READER_READER_H_
#define SRC_CONFIGURATION_READER_READER_H_

#include <string>
#include <vector>

namespace espreso {

struct GlobalConfiguration;
struct Environment;
struct Configuration;

class Reader {

public:
	static void read(Configuration &configuration, const std::string &file) { _read(configuration, file, {}); }
	static void read(Configuration &configuration, int* argc, char ***argv) { _read(configuration, argc, argv); }

	static void set(const Environment &env);

	static void print(const Configuration &configuration);
	static void store(const Configuration &configuration);

private:
	static void _read(Configuration &configuration, const std::string &file, const std::vector<std::string> &args);
	static void _read(Configuration &configuration, int* argc, char ***argv);
};

}



#endif /* SRC_CONFIGURATION_READER_READER_H_ */
