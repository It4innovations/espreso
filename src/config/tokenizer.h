
#ifndef SRC_CONFIG_TOKENIZER_H_
#define SRC_CONFIG_TOKENIZER_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "esbasis.h"

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

	Tokenizer(const std::string &file): _file(file), _token(Token::END), _line(1)
	{
		if (!_file.good()) {
			ESINFO(GLOBAL_ERROR) << "Cannot read file '" << file << "'";
		}
		_buffer.reserve(80);
	}

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



#endif /* SRC_CONFIG_TOKENIZER_H_ */
