
#ifndef SRC_CONFIGURATION_READER_TOKENIZER_H_
#define SRC_CONFIGURATION_READER_TOKENIZER_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

namespace espreso {

class Tokenizer {
public:
	enum class Token {
		STRING,
		LINK,
		DELIMITER,
		ASSIGN,
		OBJECT_OPEN,
		OBJECT_CLOSE,
		EXPRESSION_END,
		LINE_END,
		END
	};

	Tokenizer(const std::string &file);

	~Tokenizer() { _file.close(); }

	Token next() { return _token = _next(); }
	Token token() { return _token; }
	std::string value()
	{
		return std::string(_buffer.begin(), _buffer.end());
	}

	size_t line() const { return _line; }
	std::string lastLines(size_t number);

protected:
	Token _next();
	std::ifstream _file;
	Token _token;
	std::vector<int> _buffer;
	size_t _line;
};

}



#endif /* SRC_CONFIGURATION_READER_TOKENIZER_H_ */
