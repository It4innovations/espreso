
#include "tokenizer.h"

#include <queue>
#include <sstream>
#include <iomanip>

#include "../basis/logging/logging.h"

using namespace espreso;

Tokenizer::Tokenizer(const std::string &file): _file(file), _token(Token::END), _line(1)
{
	if (!_file.good()) {
		ESINFO(GLOBAL_ERROR) << "Cannot read file '" << file << "'";
	}
	_buffer.reserve(80);
}

static bool isWhiteSpace(int c)
{
	return c == ' ' || c == '\t' || c == '\r';
}

static bool isDelimiter(int c)
{
	return c == ',';
}

static bool isAssign(int c)
{
	return c == '=';
}

static bool isExpressionEnd(int c)
{
	return c == ';';
}

static bool isLineEnd(int c)
{
	return c == '\n';
}

static bool isObjectOpen(int c)
{
	return c == '{';
}

static bool isObjectClose(int c)
{
	return c == '}';
}

static bool isLinkStart(int c)
{
	return c == '[';
}

static bool isLinkEnd(int c)
{
	return c == ']';
}

static bool isStringStart(int c)
{
	return c == '"' || c == '\'';
}

static bool isStringEnd(int c)
{
	return c == '"' || c == '\'';
}

static bool isSingleCommentChar(int c)
{
	return c == '#' || c == '|';
}

static bool isLineComment(int last, int currect)
{
	return last == '/' && currect == '/';
}

static bool isMultiLineCommentStart(int last, int currect)
{
	return last == '/' && currect == '*';
}

static bool isMultiLineCommentEnd(int last, int currect)
{
	return last == '*' && currect == '/';
}

static void skipWhiteSpaces(std::ifstream &is)
{
	while (isWhiteSpace(is.peek()) && !is.eof()) {
		is.get();
	}
}

Tokenizer::Token Tokenizer::_next()
{
	if (!_file.is_open()) {
		return Token::END;
	}
	if (_file.eof()) {
		return Token::END;
	}
	skipWhiteSpaces(_file);
	_buffer.clear();

	auto specialToken = [&] (Token token) {
		if (!_buffer.size()) {
			_buffer.push_back(_file.get());
			return token;
		}
		return Token::STRING;
	};

	auto skipLine = [&] () {
		while(_file.peek() != '\n' && _file.peek() != '\r') {
			_file.get();
		}
		while(_file.peek() == '\n' || _file.peek() == '\r') {
			if (isLineEnd(_file.peek())) {
				_line++;
			}
			_file.get();
		}
	};

	while (!_file.eof()) {
		if (isSingleCommentChar(_file.peek())) {
			skipLine();
			if (_buffer.size()) {
				return Token::STRING;
			}
			skipWhiteSpaces(_file);
			continue;
		}
		if (_buffer.size() && isLineComment(_buffer.back(), _file.peek())) {
			_buffer.pop_back();
			skipLine();
			if (_buffer.size()) {
				return Token::STRING;
			}
			skipWhiteSpaces(_file);
			continue;
		}
		if (_buffer.size() && isMultiLineCommentStart(_buffer.back(), _file.peek())) {
			_buffer.pop_back();
			_file.get();
			int last = _file.get();
			while(!isMultiLineCommentEnd(last, _file.peek())) {
				if (isLineEnd(_file.peek())) {
					_line++;
				}
				last = _file.get();
			}
			_file.get();
			skipWhiteSpaces(_file);
			continue;
		}
		if (_buffer.size() && isWhiteSpace(_file.peek())) {
			_file.get();
			return Token::STRING;
		}
		if (isDelimiter(_file.peek())) {
			return specialToken(Token::DELIMITER);
		}
		if (isAssign(_file.peek())) {
			return specialToken(Token::ASSIGN);
		}
		if (isExpressionEnd(_file.peek())) {
			return specialToken(Token::EXPRESSION_END);
		}
		if (isLineEnd(_file.peek())) {
			_line++;
			return specialToken(Token::LINE_END);
		}
		if (isObjectOpen(_file.peek())) {
			return specialToken(Token::OBJECT_OPEN);
		}
		if (isObjectClose(_file.peek())) {
			return specialToken(Token::OBJECT_CLOSE);
		}
		if (isLinkStart(_file.peek())) {
			_file.get();
			while(!isLinkEnd(_file.peek())) {
				_buffer.push_back(_file.get());
			}
			_file.get();
			return Token::LINK;
		}
		if (isStringStart(_file.peek())) {
			size_t stacked = 0, bsize = _buffer.size();
			if (bsize) {
				_buffer.push_back(_file.get());
			} else {
				_file.get();
			}
			while (!isStringEnd(_file.peek()) || stacked) {
				if (_file.eof()) {
					ESINFO(GLOBAL_ERROR) << "Configuration file error: missing string end character '\"' for:\n " << std::string(_buffer.begin(), _buffer.begin() + 20);
				}
				if (isStringStart(_file.peek())) {
					stacked++;
				}
				if (isStringEnd(_file.peek())) {
					stacked--;
				}
				_buffer.push_back(_file.get());
			}
			if (bsize) {
				_buffer.push_back(_file.get());
			} else {
				_file.get();
			}
			if (_buffer.size()) {
				return Token::STRING;
			}
		}
		_buffer.push_back(_file.get());
	}

	return _buffer.size() ? Token::STRING : Token::END;
}

std::string Tokenizer::lastLines(size_t number)
{
	std::streampos current = _file.tellg();
	std::queue<std::string> lines;
	size_t line = 1;
	_file.seekg(0, _file.beg);
	while (_file.tellg() < current) {
		lines.push("");
		getline(_file, lines.back());
		if (lines.size() > number) {
			lines.pop();
		}
		line++;
	}

	line = line < number ? 0 : line - number;
	std::stringstream ss;
	while (lines.size()) {
		ss << std::setw(4) << line++ << ":" << lines.front() << "\n";
		lines.pop();
	}
	ss << std::setw(4) << line++ << ":";
	if (_file.tellg() > current) {
		return ss.str();
	}
	for (size_t i = 0; i < current - _file.tellg() + lines.size() ? lines.back().size() : 0 && i < 200; i++) {
		ss << " ";
	}
	ss << "^";
	return ss.str();
}





