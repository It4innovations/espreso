
#ifndef SRC_CONFIG_READER_READER_H_
#define SRC_CONFIG_READER_READER_H_

#include <string>
#include <vector>
#include <map>
#include <ostream>

namespace espreso {

struct ECFParameter;
struct ECFObject;
struct EnvironmentConfiguration;
struct OutputConfiguration;
struct VerboseArg;

class ECFReader {

public:
	static void read(
			ECFObject &configuration,
			const std::string &file,
			const std::map<size_t, std::string> &defaultArgs = {},
			const std::map<std::string, std::string> &variables = {})
	{
		_read(configuration, file, {}, defaultArgs, variables);
	}

	static std::string read(
			ECFObject &configuration,
			int* argc,
			char ***argv,
			const std::map<size_t, std::string> &defaultArgs = {},
			const std::map<std::string, std::string> &variables = {})
	{
		return _read(configuration, argc, argv, defaultArgs, variables);
	}

	static void store(const ECFObject &configuration, std::ostream &os, bool onlyAllowed = true, bool printPatterns = false);

private:
	static void _read(
			ECFObject &configuration,
			const std::string &file,
			const std::vector<std::string> &args,
			const std::map<size_t, std::string> &defaultArgs,
			const std::map<std::string, std::string> &variables);

	static std::string _read(
			ECFObject &configuration,
			int* argc,
			char ***argv,
			const std::map<size_t, std::string> &defaultArgs,
			const std::map<std::string, std::string> &variables);
};

}

#endif /* SRC_CONFIG_READER_READER_H_ */
