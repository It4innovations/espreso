
#include "reader.h"
#include "esbasis.h"

using namespace espreso;


void Reader::_read(Tokenizer &&tokenizer)
{
	std::vector<std::string> values;
	size_t stack = 0;

	auto indent = [] (size_t size) {
		std::stringstream _indent;
		for (size_t i = 0; i < size; i++) {
			_indent << "  ";
		}
		return _indent.str();
	};

	auto print = [&] (std::vector<std::string> &values) {
		std::cout << indent(stack);
		for (size_t i = 0; i < values.size(); i++) {
			std::cout << values[i];
		}
		std::cout << ";\n";
		values.clear();
	};

	while (tokenizer.next() != Tokenizer::Token::END) {
		switch (tokenizer.token()) {
		case Tokenizer::Token::STRING:
			values.push_back("[" + tokenizer.value() + "]");
			break;
		case Tokenizer::Token::OBJECT_OPEN:
			if (values.size() == 0) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Opening of an unnamed region is not allowed.\n" << tokenizer.lastLines(2);
			}
			if (values.size() > 1) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Multiple names for a region are not allowed.\n" << tokenizer.lastLines(2);
			}
			print(values);
			std::cout << indent(stack++) << "{\n";
			break;
		case Tokenizer::Token::OBJECT_CLOSE:
			if (!stack) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected region end.\n" << tokenizer.lastLines(2);
			}
			std::cout << indent(--stack) << "}\n";
			break;
		case Tokenizer::Token::ASSIGN:
		case Tokenizer::Token::DELIMITER:
			// assign and delimiters tokens are skipped
			break;
		case Tokenizer::Token::LINE_END:
			if (values.size() != 2) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Incorrect assignment format. Use 'PARAMETER' 'VALUE';\n" << tokenizer.lastLines(2);
			}
			print(values);
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unknown token in configuration file";
		}
	}
	if (stack) {
		ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected EOF before close all regions.";
	}

}



