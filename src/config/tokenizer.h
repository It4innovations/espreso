
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
		DELIMITER,
		ASSIGN,
		OBJECT_OPEN,
		OBJECT_CLOSE,
		EXPRESSION_END,
		LINE_END,
		END
	};

	Tokenizer(const std::string &file): _file(file), _token(Token::END)
	{
		if (!_file.good()) {
			ESINFO(GLOBAL_ERROR) << "Invalid file '" << file << "'";
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

	std::string lastLines(size_t number);

protected:
	Token _next();
	std::ifstream _file;
	Token _token;
	std::vector<int> _buffer;
};

}



#endif /* SRC_CONFIG_TOKENIZER_H_ */
