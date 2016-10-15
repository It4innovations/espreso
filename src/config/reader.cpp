
#include "reader.h"

#include <getopt.h>
#include <stack>

#include "esbasis.h"
#include "description.h"

using namespace espreso;

GlobalConfiguration espreso::configuration;

static struct option long_options[] = {
		{"configuration",  no_argument, 0, 'c'},
		{"default",  no_argument, 0, 'd'},
		{"help",  no_argument, 0, 'h'},
		{0, 0, 0, 0}
};

static std::string spaces(size_t size) {
	std::stringstream _indent;
	for (size_t i = 0; i < size; i++) {
		_indent << " ";
	}
	return _indent.str();
};

static std::string uppercase(const std::string &str) {
	std::string upper = str;
	for (auto & c: upper) { c = toupper(c); }
	return upper;
};

void Reader::_read(int* argc, char ***argv)
{
	int option_index, option;
	std::string options("c:dhvtm");

	std::vector<struct option> opts;
	std::vector<std::pair<std::string, ParameterBase*> > parameters;
	std::vector<std::string> nameless;

	std::function<void(const Configuration &conf, std::vector<std::string> path)>
	recurse = [&] (const Configuration &conf, std::vector<std::string> path) {
		for (auto it = conf.parameters.begin(); it != conf.parameters.end(); ++it) {
			std::string prefix;
			std::for_each(path.begin(), path.end(), [&] (const std::string &p) { prefix += p + "::"; });
			parameters.push_back(std::make_pair(prefix + uppercase(it->first), it->second));
		}
		for (auto it = conf.subconfigurations.begin(); it != conf.subconfigurations.end(); ++it) {
			path.push_back(uppercase(it->first));
			recurse(*it->second, path);
			path.pop_back();
		}
	};
	recurse(configuration, {});

	opts.reserve(parameters.size() + 3);
	for (size_t i = 0; i < parameters.size(); i++) {
		opts.push_back({ parameters[i].first.c_str(), required_argument, 0, 'p' });
	}

	option_index = 0;
	while (long_options[option_index].name != '\0') {
		opts.push_back(long_options[option_index++]);
	}

	// read the rest parameters
	size_t helpVerboseLevel = 0;
	std::string confFile = "espreso.ecf";
	while ((option = getopt_long(*argc, *argv, "c:dhvtm", opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'p':
			// parameters will be read after configuration file
			break;
		case 'v':
			config::info::VERBOSE_LEVEL++;
			break;
		case 't':
			config::info::TESTING_LEVEL++;
			break;
		case 'm':
			config::info::MEASURE_LEVEL++;
			break;
		case 'h':
			helpVerboseLevel++;
			break;
		case 'd':
			store();
			exit(EXIT_SUCCESS);
		case 'c':
			confFile = optarg;
			break;
		}
	}

	if (helpVerboseLevel) {
		std::cout << "PRINT HELP\n";
		exit(0);
	}

	// read nameless parameters
	while (optind < *argc) {
		nameless.push_back(std::string((*argv)[optind++]));
	}

	_read(confFile, nameless);

	optind = 0;
	while ((option = getopt_long(*argc, *argv, "c:dhvtm", opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'p':
			if (!parameters[option_index].second->set(optarg)) {
				ESINFO(GLOBAL_ERROR) << "Parameter '" << parameters[option_index].first << "' has wrong value '" << optarg << "'";
			}
			break;
		}
	}
}

void Reader::_read(const std::string &file, const std::vector<std::string> &args)
{
	std::vector<std::string> values;
	std::stack<Configuration*> confStack;
	std::stack<Tokenizer*> tokenStack;

	auto indent = [] (size_t size) {
		std::stringstream _indent;
		for (size_t i = 0; i < size; i++) {
			_indent << "  ";
		}
		return _indent.str();
	};

	confStack.push(&configuration);
	tokenStack.push(new Tokenizer(file));
	while (tokenStack.size()) {
		switch (tokenStack.top()->next()) {
		case Tokenizer::Token::END:
			delete tokenStack.top();
			tokenStack.pop();
			break;
		case Tokenizer::Token::STRING:
			if (tokenStack.top()->value()[0] == '<' && tokenStack.top()->value().back() == '>') {
				std::string value = tokenStack.top()->value();
				value = value.substr(1, value.size() - 2);
				if (value.size() > 4 && StringCompare::caseInsensitiveSuffix(value, ".ecf")) {
					tokenStack.push(new Tokenizer(value));
				}
				if (value.size() > 2 && StringCompare::caseInsensitivePreffix("ARG", value)) {
					std::stringstream ss(std::string(value.begin() + 3, value.end()));
					int index = args.size();
					ss >> index;
					if (!ss.fail() && ss.eof() && index < args.size()) {
						values.push_back(args[index]);
					} else {
						if (index < args.size()) {
							ESINFO(GLOBAL_ERROR) << "Invalid argument '" << value << "'";
						} else {
							ESINFO(GLOBAL_ERROR) << "ARG index is out of range '" << value << "'";
						}
					}
				}
			} else {
				values.push_back(tokenStack.top()->value());
			}
			break;
		case Tokenizer::Token::OBJECT_OPEN:
			if (values.size() == 0) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Opening of an unnamed region is not allowed.\n" << tokenStack.top()->lastLines(2);
			}
			if (values.size() > 1) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Multiple names for a region are not allowed.\n" << tokenStack.top()->lastLines(2);
			}
			confStack.push(&confStack.top()->operator [](values[0]));
			values.clear();
			break;
		case Tokenizer::Token::OBJECT_CLOSE:
			if (!confStack.size()) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected region end.\n" << tokenStack.top()->lastLines(2);
			}
			confStack.pop();
			break;
		case Tokenizer::Token::ASSIGN:
		case Tokenizer::Token::DELIMITER:
			// assign and delimiters tokens are skipped
			break;
		case Tokenizer::Token::EXPRESSION_END:
			if (values.size() != 2) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Incorrect assignment format. Use 'PARAMETER' 'VALUE';\n" << tokenStack.top()->lastLines(2);
			}
			if (!confStack.top()->set(values[0], values[1])) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Parameter '" << values[0] << "' has wrong value '" << values[1] << "'";
			}
			values.clear();
			break;
		case Tokenizer::Token::LINE_END:
			if (values.size() > 1) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Expected ';' at the end of each expression.\n" << tokenStack.top()->lastLines(1);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unknown token in configuration file";
		}
	}
	if (confStack.size() != 1) {
		ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected EOF before close all regions.";
	}

}

static void printConfiguration(const Configuration &configuration, size_t indent)
{
	for (auto it = configuration.parameters.begin(); it != configuration.parameters.end(); ++it) {
		ESINFO(ALWAYS) << spaces(indent) << uppercase(it->first) << " = " << uppercase(it->second->get());
	}

	for (auto it = configuration.subconfigurations.begin(); it != configuration.subconfigurations.end(); ++it) {
		ESINFO(ALWAYS) << spaces(indent) << uppercase(it->first) << " {";
		printConfiguration(*it->second, indent + 2);
		ESINFO(ALWAYS) << spaces(indent) << "}";
	}
}

static void storeConfiguration(std::ofstream &os, const Configuration &configuration, size_t indent)
{
	for (auto it = configuration.parameters.begin(); it != configuration.parameters.end(); ++it) {
		os << "\n" << spaces(indent) << "# " << it->second->description << "\n";
		os << spaces(indent) << uppercase(it->first) << " " << uppercase(it->second->get()) << ";\n";
	}

	for (auto it = configuration.subconfigurations.begin(); it != configuration.subconfigurations.end(); ++it) {
		os << "\n" << spaces(indent) << uppercase(it->first) << " {\n";
		storeConfiguration(os, *it->second, indent + 2);
		os << spaces(indent) << "}\n\n";
	}
}

void Reader::print()
{
	ESINFO(ALWAYS) << "ESPRESO configuration:";
	printConfiguration(configuration, 4);
}

void Reader::store()
{
	std::ofstream os("espreso.ecf.default");

	os << "*******************************************************************************\n";
	os << "|-----------------------------------------------------------------------------|\n";
	os << "|                                      |                                      |\n";
	os << "|     EPRESO CONFIGURATION FILE        |   ESPRESO Version:   1.0             |\n";
	os << "|                                      |   http://espreso.it4i.cz             |\n";
	os << "|-----------------------------------------------------------------------------|\n";
	os << "|  Case Description:    Default ESPRESO configuration                         |\n";
	os << "|                                                                             |\n";
	os << "|-----------------------------------------------------------------------------|\n";
	os << "*******************************************************************************\n";
	os << "                                                                               \n";
	os << "                                                                               \n";
	os << "*******************************************************************************\n";
	os << "|-------------------------  INPUT/OUTPUT DEFINITION --------------------------|\n\n";

	storeConfiguration(os, configuration, 0);
	ESINFO(ALWAYS) << "configuration stored to 'espreso.ecf.default'";
}



