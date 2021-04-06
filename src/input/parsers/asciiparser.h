
#ifndef SRC_INPUT_PARSERS_ASCIIPARSER_H_
#define SRC_INPUT_PARSERS_ASCIIPARSER_H_

#include <vector>
#include <cstddef>

namespace espreso {

struct InputFile;

class ASCIIParser {
public:
	static inline bool isempty(const char *c)
	{
		return *c == ' ' || *c == '\n' || *c == '\t' || *c == '\r';
	}
	static int keyend(const char *c);

	static void parse(std::vector<esint> &data, InputFile &file, size_t begin, size_t end);
	static void parse(std::vector<double> &data, InputFile &file, size_t begin, size_t end);

	static void addmore(std::vector<esint> &data, InputFile &file, size_t n, size_t end);
	static void addmore(std::vector<double> &data, InputFile &file, size_t n, size_t end);

private:
	template <typename TType>
	static void _parse(std::vector<TType> &data, InputFile &file, size_t begin, size_t end);
	template <typename TType>
	static void _addmore(std::vector<TType> &data, InputFile &file, size_t n, size_t end);
};

}

#endif /* SRC_INPUT_PARSERS_ASCIIPARSER_H_ */
