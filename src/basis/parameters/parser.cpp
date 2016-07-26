
#include "parser.h"

using namespace espreso;

std::string Parser::strip(const std::string &line)
{
	std::string result(line);
	result.erase(0, line.find_first_not_of(" \r\t\n"));
	result.erase(result.find_last_not_of(" \r\t\n") + 1, std::string::npos);
	return result.substr(0, result.find("#"));
}

std::vector<std::string> Parser::split(const std::string &line, const std::string &separator)
{
	if (line.size() == 0) {
		return std::vector<std::string> (1);
	}

	std::vector<std::string> result;
	std::string reminder = line;

	while (reminder.size()) {
		result.push_back(reminder.substr(0, reminder.find_first_of(separator)));
		if (reminder.find_first_of(separator) == std::string::npos) {
			reminder.erase(0, std::string::npos);
		} else {
			reminder.erase(0, reminder.find_first_of(separator) + 1);
		}
	}
	return result;
}

bool Parser::contains(const std::string &line, const std::string &separator)
{
	return line.find(separator) != std::string::npos;
}

std::string Parser::getParameter(const std::string &line, const std::string &separator)
{
	return strip(split(strip(line), separator)[0]);
}

std::string Parser::getValue(const std::string &line, const std::string &separator)
{
	if (split(strip(line), separator).size() > 1) {
		return strip(split(strip(line), separator)[1]);
	}
	return "";
}




