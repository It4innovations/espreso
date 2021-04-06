
#ifndef SRC_CONFIG_READER_TOKENIZER_H_
#define SRC_CONFIG_READER_TOKENIZER_H_

#include <string>
#include <vector>

namespace espreso {

class FullFile {

public:
	FullFile(const std::string &file);

	size_t position() const { return _p; }
	void move(size_t position) { _p = position; }

	char get() { return _data[_p++]; }
	char peek() { return _data[_p]; }

	std::string get(size_t begin, size_t end) { return std::string(_data.begin() + begin, _data.begin() + end); }

	bool eof() const { return _p == _data.size(); }

	void synchronize();

protected:
	std::vector<char> _data;
	size_t _p;
};

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

	FullFile _file;
	Token _token;
	std::vector<int> _buffer;
	size_t _line;
};

class CollectiveTokenizer: public Tokenizer {
public:
	CollectiveTokenizer(const std::string &file): Tokenizer(file)
	{
		_file.synchronize();
	}
};

}

#endif /* SRC_CONFIG_READER_TOKENIZER_H_ */
