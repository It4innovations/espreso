
#ifndef SRC_BASIS_PARAMETERS_READER_H_
#define SRC_BASIS_PARAMETERS_READER_H_

#include "mpi.h"

#include <getopt.h>

#include "esconfig.h"
#include "parameter.h"
#include "../logging/logging.h"


namespace espreso {

class ParametersReader {

public:
	static Configuration fromArguments(int *argc, char*** argv, const std::vector<Parameter> &params = config::parameters);
	static Configuration fromConfigurationFile(const Configuration &conf, const std::vector<Parameter> &params = config::parameters);
	static Configuration fromConfigurationFileWOcheck(const Configuration &conf, const std::vector<Parameter> &params = config::parameters);

	static void printParameters(const std::vector<Parameter> &params, size_t verboseLevel);
	static void printParametersHelp(const std::vector<Parameter> &params, size_t verboseLevel);

	ParametersReader(const std::vector<Parameter> &parameters);

	template<typename Tvalue>
	bool setParameter(const void* parameter, Tvalue value)
	{
		std::stringstream ss;
		ss << value;
		for (size_t i = 0; i < _parameters.size(); i++) {
			if (_parameters[i].data->value() == parameter) {
				return _parameters[i].data->set(ss.str());
			}
		}
		return false;
	}

protected:
	Configuration read(const Configuration &configuration, size_t verboseLevel);

	std::vector<Parameter> _parameters;

private:
	static void printHelp(size_t verboseLevel);
	bool _setParameter(const std::string &parameter, const std::string &value);
	bool _setParameters(const std::string &parameter, const std::vector<std::string> &attributes, const std::vector<std::string> &values);
};

}



#endif /* SRC_BASIS_PARAMETERS_READER_H_ */
