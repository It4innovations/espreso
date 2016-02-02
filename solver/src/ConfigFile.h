#include "utils.h"

#include <string>
#include <map>
#include <vector>



class ConfigFile {
	std::map   <std::string, std::string> content_;
	std::vector<std::string>              sections_;
public:
	ConfigFile(std::string const& configFile);

	std::vector<std::string> GetSections();
	std::string const& Value(std::string const& section, std::string const& entry) const;
	std::string const& Value(std::string const& section, std::string const& entry, std::string const& value);
};
