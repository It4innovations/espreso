
#include "parser.h"

using namespace espreso;

std::string Parser::getLine(const char* begin)
{
	const char* end = begin;
	while (*end++ != '\n');
	return std::string(begin, end);
}

std::string Parser::uppercase(const std::string &str) {
	std::string upper = str;
	for (auto & c: upper) { c = toupper(c); }
	return upper;
};

std::string Parser::removecomments(const std::string &line, const std::string &tags)
{
	return line.substr(0, line.find_first_of(tags));
}

std::string Parser::strip(const std::string &line)
{
	std::string result(line);
	result.erase(0, line.find_first_not_of(" \r\t\n"));
	result.erase(result.find_last_not_of(" \r\t\n") + 1, std::string::npos);
	return result.substr(0, result.find("#"));
}

std::vector<std::string> Parser::split(const std::string &line, const std::string &separator, bool skipMultiple)
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
			if (skipMultiple) {
				reminder.erase(0, reminder.find_first_of(separator));
				reminder.erase(0, reminder.find_first_not_of(separator));
			} else {
				reminder.erase(0, reminder.find_first_of(separator) + 1);
			}
		}
	}
	return result;
}

std::vector<std::pair<std::string, std::string> > Parser::getIntervals(const std::string &line)
{
	if (line.size() == 0) {
		return std::vector<std::pair<std::string, std::string> >();
	}

	std::vector<std::pair<std::string, std::string> > result;

	size_t pos = 0;
	while (pos < line.size()) {
		size_t begin = line.find_first_of('<', pos);
		size_t mid = line.find_first_of(',', pos);
		size_t end = line.find_first_of('>', pos);
		if (begin < line.size() && end < line.size() && mid < line.size() && begin < mid && mid < end) {
			result.push_back(std::make_pair(line.substr(begin + 1, mid - begin - 1), line.substr(mid + 1, end - mid - 1)));
			pos = end + 1;
		} else {
			pos = line.size();
		}
	}
	return result;
}

std::string Parser::stringwithcommas(size_t number)
{
	std::string result;
	long value = 1e9;
	bool zeros = false;
	for (int i = 0; i < 3; ++i) {
		if (zeros) {
			if (number / value) {
				result += std::string(3 - std::to_string(number / value).size(), '0');
			} else {
				result += "000,";
			}
		}
		if (number / value) {
			result += std::to_string(number / value) + ",";
			zeros = true;
		}
		number %= value;
		value /= 1000;
	}
	if (zeros) {
		result += std::string(3 - std::to_string(number).size(), '0');
	}
	result += std::to_string(number);
	return result;
};



