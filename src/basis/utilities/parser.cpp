
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

std::string Parser::join(const std::string &separator, const std::vector<std::string> &values)
{
    if (values.size() == 0) {
        return "";
    }
    std::string result = values.front();
    for (size_t i = 1; i < values.size(); ++i) {
        result += separator + values[i];
    }
    return result;
}



