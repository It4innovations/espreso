
#ifndef SRC_BASIS_PARAMETERS_PARSER_H_
#define SRC_BASIS_PARAMETERS_PARSER_H_

#include <string>
#include <vector>
#include <iostream>

namespace espreso {

struct StringCompare {

	bool operator()(const std::string &s1, const std::string &s2)
	{
		return caseInsensitive(s1, s2);
	}

	static bool caseSensitiveEq(const std::string &s1, const std::string &s2)
	{
		return s1.size() == s2.size() && std::equal(s1.begin(), s1.end(), s2.begin());
	}

	static bool caseInsensitiveEq(const std::string &s1, const std::string &s2)
	{
		return s1.size() == s2.size() && std::equal(s1.begin(), s1.end(), s2.begin(), _equals);
	}

	static bool caseSensitive(const std::string &s1, const std::string &s2)
	{
		return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end());
	}

	static bool caseInsensitive(const std::string &s1, const std::string &s2)
	{
		return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end(), _caseInsensitive);
	}

	static bool caseSensitivePreffix(const std::string &s1, const std::string &s2)
	{
		return std::equal(s1.begin(), s1.end(), s2.begin());
	}

	static bool caseInsensitivePreffix(const std::string &s1, const std::string &s2)
	{
		return std::equal(s1.begin(), s1.end(), s2.begin(), _equals);
	}

	static bool caseSensitiveSuffix(const std::string &s1, const std::string &s2)
	{
		return std::equal(s1.begin() + s1.size() - s2.size(), s1.end(), s2.begin());
	}

	static bool caseInsensitiveSuffix(const std::string &s1, const std::string &s2)
	{
		return std::equal(s1.begin() + s1.size() - s2.size(), s1.end(), s2.begin(), _equals);
	}

private:
	static bool _caseInsensitive(const char &c1, const char &c2)
	{
		return std::tolower(c1) < std::tolower(c2);
	}

	static bool _equals(const char &c1, const char &c2)
	{
		return std::tolower(c1) == std::tolower(c2);
	}
};

class Parser {
public:
	static std::string getParameter(const std::string &line, const std::string &separator = "=");
	static std::string getValue(const std::string &line, const std::string &separator = "=");
	static bool contains(const std::string &line, const std::string &separator);

	static std::string strip(const std::string &line);
	static std::vector<std::string> split(const std::string &line, const std::string &separator);
};

}



#endif /* SRC_BASIS_PARAMETERS_PARSER_H_ */
