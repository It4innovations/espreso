
#include "utils.h"

#include <fstream>
#include <algorithm>

using namespace espreso;

std::string Esutils::createDirectory(const std::vector<std::string> &path)
{
	std::stringstream prefix;
	std::for_each(path.begin(), path.end(), [&] (const std::string &dir) { prefix << dir << "/"; });

	if (system(("mkdir -p " + prefix.str()).c_str())) {
		ESINFO(ERROR) << "Cannot create requested directory";
	}
	return prefix.str();
}


