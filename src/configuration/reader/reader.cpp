
#include "reader.h"

#include <getopt.h>
#include <stack>
#include <functional>
#include <regex>

#include "../basis/logging/logging.h"
#include "../configuration/globalconfiguration.h"
#include "tokenizer.h"

using namespace espreso;

static struct option long_options[] = {
		{"config",  required_argument, 0, 'c'},
		{"default",  optional_argument, 0, 'd'},
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

void Reader::_read(Configuration &configuration, int* argc, char ***argv)
{
	environment->executable = *argv[0];
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
	if (StringCompare::caseSensitiveSuffix(std::string(*argv[0]), "espreso")) {
		confFile = "espreso.ecf";
	}
	if (StringCompare::caseSensitiveSuffix(std::string(*argv[0]), "decomposer")) {
		confFile = "decomposer.ecf";
	}

	std::vector<std::string> subConfigurations;
	while ((option = getopt_long(*argc, *argv, "c:d::hvtm", opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'p':
			// parameters will be read after configuration file
			break;
		case 'h':
			helpVerboseLevel++;
			break;
		case 'd':
			if (optarg == NULL) {
				subConfigurations.push_back(".*");
			} else {
				subConfigurations.push_back(std::string(optarg));
			}
			break;
		case 'c':
			confFile = optarg;
			break;
		case '?':
			exit(EXIT_FAILURE);
			break;
		}
	}

	if (subConfigurations.size()) {
		store(configuration, subConfigurations);
		exit(EXIT_SUCCESS);
	}

	if (helpVerboseLevel) {
		std::cout << "\nusage: ./espreso [options] [ARGS]\n\n";

		std::cout << " [options] are the following:\n";
		std::cout << "\t -h, --help            print this message\n";
		std::cout << "\t -d, --default         generate default configuration file\n";
		std::cout << "\t -c, --config=[path]   set path to configuration file\n\n";

		std::cout << " [ARGS] are unnamed argument that can be referenced in a configuration\n";
		std::cout << "        file by [ARGX], where X is a number of an argument counted from 0.\n";
		std::cout << "        e.g. DOMAINS [ARG0]; set number of domains to the first argument.\n\n";

		std::cout << "The solver is controlled by '*.ecf' scripts. Some examples can be found\n";
		std::cout << "in 'benchmarks' directory or in 'tests/examples'. In default 'espreso.ecf'\n";
		std::cout << "from ESPRESO root directory is loaded. A different configuration file can\n";
		std::cout << "be set by -c [path] option.\n\n";
		std::cout << "A file is composed from parameters and objects with the following pattern:\n\n";
		std::cout << "  OBJECT {\n";
		std::cout << "    PARAMETER VALUE;\n";
		std::cout << "  }\n\n";
		std::cout << "Run ./espreso --default to generates the configuration file with default\n";
		std::cout << "parameters.\n\n";
		exit(0);
	}

	// read nameless parameters
	while (optind < *argc) {
		nameless.push_back(std::string((*argv)[optind++]));
	}

	_read(configuration, confFile, nameless);

	optind = 0;
	while ((option = getopt_long(*argc, *argv, "c:dhvtm", opts.data(), &option_index)) != -1) {
		switch (option) {
		case 'p':
			if (!parameters[option_index].second->set(optarg)) {
				ESINFO(GLOBAL_ERROR) << "Parameter '" << parameters[option_index].first << "' has wrong value '" << optarg << "'";
			}
			break;
		case 'v':
			environment->verbose_level++;
			break;
		case 't':
			environment->testing_level++;
			break;
		case 'm':
			environment->measure_level++;
			break;
		}
	}
}

void Reader::_read(Configuration &configuration, const std::string &file, const std::vector<std::string> &args)
{
	std::vector<std::string> prefix;
	std::vector<std::string> values;
	std::stack<Configuration*> confStack;
	std::stack<Tokenizer*> tokenStack;

	bool correctlyLoaded = true;
	std::map<size_t, std::vector<std::string> > arguments;

	confStack.push(&configuration);
	tokenStack.push(new Tokenizer(file));
	while (tokenStack.size()) {
		switch (tokenStack.top()->next()) {
		case Tokenizer::Token::END:
			delete tokenStack.top();
			tokenStack.pop();
			break;
		case Tokenizer::Token::STRING:
			values.push_back(tokenStack.top()->value());
			break;
		case Tokenizer::Token::LINK:
		{
			std::string value = tokenStack.top()->value();
			if (value.size() > 4 && StringCompare::caseInsensitiveSuffix(value, ".ecf")) {
				tokenStack.push(new Tokenizer(value));
				break;
			}
			if (value.size() > 2 && StringCompare::caseInsensitivePreffix("ARG", value)) {
				std::stringstream ss(std::string(value.begin() + 3, value.end()));
				size_t index = args.size();
				ss >> index;
				if (!ss.fail() && ss.eof() && index < args.size()) {
					values.push_back(args[index]);
					std::string parameter;
					std::for_each(prefix.begin(), prefix.end(), [&] (const std::string &s) { parameter += s + "::"; });
					arguments[index].push_back(parameter + values.front());
				} else {
					if (index < args.size()) {
						ESINFO(GLOBAL_ERROR) << "Invalid argument '" << value << "'";
					} else {
						correctlyLoaded = false;
						if (values.size()) {
							std::string parameter;
							std::for_each(prefix.begin(), prefix.end(), [&] (const std::string &s) { parameter += s + "::"; });
							arguments[index].push_back(parameter + values.front());
						} else {
							ESINFO(GLOBAL_ERROR) << "parameter cannot be the [ARG].\n" << tokenStack.top()->lastLines(2);
						}
					}
				}
				break;
			}
			if (value.size() > 4 && StringCompare::caseInsensitiveSuffix(value, ".csv")) {
				Tokenizer csv(value);
				bool run = true;
				while(run) {
					switch (csv.next()) {
					case Tokenizer::Token::END:
						run = false;
						break;
					case Tokenizer::Token::STRING:
						values.push_back(csv.value() + (values.size() % 2 == 0 ? "," : ";"));
						break;
					case Tokenizer::Token::LINE_END:
					case Tokenizer::Token::DELIMITER:
					case Tokenizer::Token::EXPRESSION_END:
						break;
					default:
						ESINFO(GLOBAL_ERROR) << "Error while reading file '" << value << "'";
					}
				}
				break;
			}
			values.push_back(value);
			break;
		}
		case Tokenizer::Token::OBJECT_OPEN:
			if (values.size() == 0) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Opening of an unnamed region is not allowed.\n" << tokenStack.top()->lastLines(2);
			}
			if (values.size() > 1) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Multiple names for a region are not allowed.\n" << tokenStack.top()->lastLines(2);
			}
			prefix.push_back(values[0]);
			confStack.push(&confStack.top()->operator [](values[0]));
			values.clear();
			break;
		case Tokenizer::Token::OBJECT_CLOSE:
			if (!confStack.size()) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Unexpected region end.\n" << tokenStack.top()->lastLines(2);
			}
			prefix.pop_back();
			confStack.pop();
			break;
		case Tokenizer::Token::ASSIGN:
			break;
		case Tokenizer::Token::DELIMITER:
			values.push_back(",");
			break;
		case Tokenizer::Token::EXPRESSION_END:
		{
			if (!correctlyLoaded) {
				values.clear();
				break;
			}
			if (values.size() < 2) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Incorrect assignment format on line " << tokenStack.top()->line() << ". Use 'PARAMETER' 'VALUE';\n" << tokenStack.top()->lastLines(2);
			}
			std::stringstream ss;
			ss << values[1];
			for (size_t i = 2; i < values.size(); i++) {
				ss << " " << values[i];
			}
			if (!confStack.top()->set(values[0], ss.str())) {
				ESINFO(GLOBAL_ERROR) << "PARSE ERROR: Parameter '" << values[0] << "' has wrong value '" << ss.str() << "'";
			}

			values.clear();
			break;
		}
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

	if (!correctlyLoaded) {
		std::string error = "Configuration file is not correctly loaded.\nUse ./espreso ";
		size_t i = 0;
		for (auto it = arguments.begin(); it != arguments.end(); ++it) {
			error += "[ARG" + std::to_string(i++) + "] ";
		}
		error += "\nWhere ARGs are the following:\n";
		i = 0;
		for (auto it = arguments.begin(); it != arguments.end(); ++it) {
			error += "ARG" + std::to_string(i++) + " = { ";
			for (size_t j = 0; j < it->second.size(); j++) {
				error += it->second[j];
				if (j != it->second.size() - 1) {
					error += ", ";
				}
			}
			error += " }\n";
		}
		ESINFO(GLOBAL_ERROR) << error;
	}
}

void Reader::set(const Environment &env)
{
	Test::setLevel(env.testing_level);
	Info::setLevel(env.verbose_level, env.testing_level);
	Measure::setLevel(env.measure_level);
	Logging::output = env.log_dir;
	Logging::rank = env.MPIrank;
}

static void printConfiguration(const Configuration &configuration, size_t indent)
{
	for (size_t i = 0; i < configuration.orderedParameters.size(); i++) {
		ParameterBase *parameter = configuration.orderedParameters[i];
		ESINFO(ALWAYS) << spaces(indent) << uppercase(parameter->name) << " = " << parameter->get();
	}

	for (size_t i = 0; i < configuration.orderedSubconfiguration.size(); i++) {
		ESINFO(ALWAYS) << spaces(indent) << uppercase(configuration.orderedSubconfiguration[i]->name) << " {";
		printConfiguration(*configuration.orderedSubconfiguration[i], indent + 2);
		ESINFO(ALWAYS) << spaces(indent) << "}";
	}
}

static void storeConfiguration(std::ofstream &os, const Configuration &configuration, size_t indent, std::vector<std::regex> &patterns)
{
	for (size_t i = 0; i < configuration.storeParameters().size(); i++) {
		ParameterBase *parameter = configuration.storeParameters()[i];
		os << "\n" << spaces(indent) << "# " << parameter->description << " [" << parameter->allowedValue << "]\n";
		os << spaces(indent) << uppercase(parameter->name) << " " << parameter->get() << ";\n";
	}

	for (size_t i = 0; i < configuration.storeConfigurations().size(); i++) {
		if (std::any_of(patterns.begin(), patterns.end(), [&] (const std::regex &regex) {
			std::smatch sm;
			std::regex_match(configuration.storeConfigurations()[i]->name, sm, regex);
			return sm.size();
		})) {

			os << "\n" << spaces(indent) << uppercase(configuration.storeConfigurations()[i]->name) << " { ";
			os << "# " << configuration.storeConfigurations()[i]->description << "\n";
			std::vector<std::regex> all = { std::regex(".*") };
			storeConfiguration(os, *configuration.storeConfigurations()[i], indent + 2, all);
			os << spaces(indent) << "}\n\n";
		}
	}
}

void Reader::print(const Configuration &configuration)
{
	ESINFO(ALWAYS) << "ESPRESO configuration:";
	printConfiguration(configuration, 4);
}

void Reader::store(const Configuration &configuration, const std::vector<std::string> &subConfigurations)
{
	std::ofstream os("espreso.ecf.default");

	os << "|*****************************************************************************|\n";
	os << "|-----------------------------------------------------------------------------|\n";
	os << "|                                      |                                      |\n";
	os << "|     ESPRESO CONFIGURATION FILE       |   ESPRESO Version:   1.0             |\n";
	os << "|                                      |   http://espreso.it4i.cz             |\n";
	os << "|-----------------------------------------------------------------------------|\n";
	os << "|  Case Description:    Default ESPRESO configuration                         |\n";
	os << "|                                                                             |\n";
	os << "|-----------------------------------------------------------------------------|\n";
	os << "|*****************************************************************************|\n";
	os << "                                                                               \n";
	os << "                                                                               \n";
	os << "|*****************************************************************************|\n";
	os << "|-------------------------  INPUT/OUTPUT DEFINITION --------------------------|\n\n";

	std::vector<std::regex> patterns;
	std::for_each(subConfigurations.begin(), subConfigurations.end(), [&] (const std::string &s) { patterns.push_back(std::regex(s)); });

	storeConfiguration(os, configuration, 0, patterns);
	ESINFO(ALWAYS) << "configuration stored to 'espreso.ecf.default'";
}




