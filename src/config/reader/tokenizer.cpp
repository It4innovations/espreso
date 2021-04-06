#include "tokenizer.h"

#include "mpi.h"

#include <queue>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "basis/utilities/communication.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

using namespace espreso;

FullFile::FullFile(const std::string &file): _p(0)
{
	if (info::mpi::rank == 0) {
		std::ifstream data(file, std::ifstream::in);
		if (!data.good()) {
			eslog::error("Cannot read file '%s'\n", file.c_str());
		}

		data.seekg(0, data.end);
		_data.resize(data.tellg());
		data.seekg(0, data.beg);

		data.read(_data.data(), _data.size());
		data.close();
	}
}

void FullFile::synchronize()
{
	size_t filesize = _data.size();
	Communication::broadcast(&filesize, sizeof(size_t), MPI_BYTE, 0, MPITools::global);
	_data.resize(filesize);
	Communication::broadcast(_data.data(), _data.size(), MPI_CHAR, 0, MPITools::global);
}


Tokenizer::Tokenizer(const std::string &file): _file(file), _token(Token::END), _line(1)
{
	_buffer.reserve(80);
}

static bool isWhiteSpace(char c)
{
	return c == ' ' || c == '\t' || c == '\r';
}

static bool isDelimiter(char c)
{
	return c == ',';
}

static bool isAssign(char c)
{
	return c == '=';
}

static bool isExpressionEnd(char c)
{
	return c == ';';
}

static bool isLineEnd(char c)
{
	return c == '\n';
}

static bool isObjectOpen(char c)
{
	return c == '{';
}

static bool isObjectClose(char c)
{
	return c == '}';
}

static bool isLinkStart(char c)
{
	return c == '[';
}

static bool isLinkEnd(char c)
{
	return c == ']';
}

static bool isStringStart(char c)
{
	return c == '"' || c == '\'';
}

static bool isStringEnd(char c)
{
	return c == '"' || c == '\'';
}

static bool isSingleCommentChar(char c)
{
	return c == '#';
}

static bool isLineComment(char last, char currect)
{
	return last == '/' && currect == '/';
}

static bool isMultiLineCommentStart(char last, char currect)
{
	return last == '/' && currect == '*';
}

static bool isMultiLineCommentEnd(char last, char currect)
{
	return last == '*' && currect == '/';
}

static void skipWhiteSpaces(FullFile &file)
{
	while (isWhiteSpace(file.peek()) && !file.eof()) {
		file.get();
	}
}

Tokenizer::Token Tokenizer::_next()
{
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
		while (!_file.eof() && _file.peek() != '\n' && _file.peek() != '\r') {
			_file.get();
		}
		while (!_file.eof() && (_file.peek() == '\n' || _file.peek() == '\r')) {
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
		if (isLinkStart(_file.peek()) && _buffer.size()) {
			return Token::STRING;
		}
		if (isLinkStart(_file.peek())) {
			_file.get();
			size_t stacked = 0;
			while(!isLinkEnd(_file.peek()) || stacked) {
				if (_file.eof()) {
					eslog::globalerror("Configuration file error: missing link end symbol ']'");
				}
				if (isLinkStart(_file.peek())) {
					stacked++;
				}
				if (isLinkEnd(_file.peek())) {
					stacked--;
				}
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
					eslog::globalerror("Configuration file error: missing link end symbol '\"'");
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
			return Token::STRING;
		}
		_buffer.push_back(_file.get());
	}

	return _buffer.size() ? Token::STRING : Token::END;
}

std::string Tokenizer::lastLines(size_t number)
{
	size_t current = _file.position();
	std::queue<std::string> lines;
	size_t line = 1;
	_file.move(0);
	while (_file.position() < current) {
		size_t begin = _file.position();
		while(_file.get() != '\n');
		lines.push(_file.get(begin, _file.position()));
		if (lines.size() > number) {
			lines.pop();
		}
		line++;
	}

	line = line < number ? 0 : line - number;
	std::stringstream ss;
	while (lines.size()) {
		ss << std::setw(4) << line++ << ":" << lines.front();
		lines.pop();
	}
	ss << std::setw(4) << line++ << ":";
	if (_file.position() > current) {
		return ss.str();
	}
	for (size_t i = 0; i < (current - _file.position()) + lines.size() ? lines.back().size() : 0 && i < 200; i++) {
		ss << " ";
	}
	ss << "^";
	return ss.str();
}



