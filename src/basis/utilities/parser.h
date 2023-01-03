
#ifndef SRC_BASIS_UTILITIES_PARSER_H_
#define SRC_BASIS_UTILITIES_PARSER_H_

#include <string>
#include <vector>

namespace espreso {

struct StringCompare {

	bool operator()(const std::string &s1, const std::string &s2) const
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
		if (s1.size() <= s2.size()) {
			return std::equal(s1.begin(), s1.end(), s2.begin(), _equals);
		}
		return false;
	}

	static bool caseSensitiveSuffix(const std::string &s1, const std::string &s2)
	{
		if (s1.size() >= s2.size()) {
			return std::equal(s1.begin() + s1.size() - s2.size(), s1.end(), s2.begin());
		}
		return false;
	}

	static bool caseInsensitiveSuffix(const std::string &s1, const std::string &s2)
	{
		if (s1.size() >= s2.size()) {
			return std::equal(s1.begin() + s1.size() - s2.size(), s1.end(), s2.begin(), _equals);
		}
		return false;
	}

	static bool contains(const std::string &s, const std::vector<std::string> &variables) {
		std::string lower = s;
		for (size_t i = 0; i < lower.size(); i++) {
			lower[i] = std::tolower(lower[i]);
		}
		for (size_t v = 0; v < variables.size(); v++) {
			std::string lvariable = variables[v];
			for (size_t i = 0; i < lvariable.size(); i++) {
				lvariable[i] = std::tolower(lvariable[i]);
			}
			if (lower.find(lvariable) != std::string::npos) {
				return true;
			}
		}
		return false;
	};

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
	static std::string getLine(const char* begin);

	static std::string uppercase(const std::string &str);

	static std::string removecomments(const std::string &line, const std::string &tags);
	static std::string strip(const std::string &line);
	static std::vector<std::string> split(const std::string &line, const std::string &separator, bool skipMultiple = true);
	static std::vector<std::pair<std::string, std::string> > getIntervals(const std::string &line);

	static std::string stringwithcommas(size_t number);
};

}



#endif /* SRC_BASIS_UTILITIES_PARSER_H_ */
