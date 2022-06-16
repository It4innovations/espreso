
#include "configuration.h"
#include "conditions.h"
#include "esinfo/eslog.hpp"
#include "basis/utilities/parser.h"
#include <algorithm>

using namespace espreso;

void ECFParameter::defaultName()
{
	if (!metadata.name.size()) {
		metadata.name = name;

		for (size_t i = 0, up = 1; i < metadata.name.size(); i++) {
			if (up && 'a' <= metadata.name[i] && metadata.name[i] <= 'z') {
				metadata.name[i] += 'A' - 'a';
			}
			up = 0;
			if (metadata.name[i] == '_') {
				metadata.name[i] = ' ';
				up = 1;
			}
		}
	}
}

bool ECFParameter::isvisible()
{
	return metadata.condition->evaluate();
}


bool ECFParameter::setValue(const std::string &value)
{
	if (_setValue(value)) {
		for (size_t i = 0; i < _setValueListeners.size(); i++) {
			_setValueListeners[i](value);
		}
		return true;
	}
	return false;
}

ECFParameter* ECFParameter::_triggerParameterGet(ECFParameter* parameter)
{
	if (parameter != NULL) {
		for (size_t i = 0; i < _parameterGetListeners.size(); i++) {
			_parameterGetListeners[i](parameter->name);
		}
	}
	return parameter;
}

ECFParameter* ECFParameter::getParameter(const std::string &name)
{
	return _triggerParameterGet(_getNamedParameter(name));
}

ECFParameter* ECFParameter::getParameter(const char* name)
{
	return _triggerParameterGet(_getNamedParameter(std::string(name)));
}

ECFParameter* ECFParameter::getParameter(const void* data)
{
	return _triggerParameterGet(_getParameter(data));
}

void ECFParameter::addListener(Event event, std::function<void(const std::string &value)> listener)
{
	switch (event) {
	case Event::VALUE_SET:
		_setValueListeners.push_back(listener);
		break;
	case Event::PARAMETER_GET:
		_parameterGetListeners.push_back(listener);
	}
}

ECFParameter* ECFParameter::registerAdditionalParameter(ECFParameter* parameter)
{
	eslog::internalFailure("parameter '%s' cannot have additional parameters.\n", name.c_str());
	return NULL;
}

std::string ECFObject::getValue() const
{
	eslog::internalFailure("calling 'get' on ECFObject.\n");
	return "";
}

bool ECFObject::_setValue(const std::string &value)
{
	eslog::internalFailure("calling 'set' on ECFObject.\n");
	return false;
}

ECFParameter* ECFObject::_getNamedParameter(const std::string &name)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (StringCompare::caseInsensitiveEq(name, parameters[i]->name)) {
			return parameters[i];
		}
	}
	return NULL;
}

ECFParameter* ECFObject::_getParameter(const void* data)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (data == parameters[i]->data()) {
			return parameters[i];
		}
	}
	return NULL;
}

ECFParameter* ECFObject::getWithError(const std::string &name)
{
	if (_getNamedParameter(name) == NULL) {
		eslog::globalerror("ECF ERROR: Object '%s' has no parameter '%s'.\n", this->name.c_str(), name.c_str());
	}
	return _getNamedParameter(name);
}

ECFParameter* ECFObject::registerAdditionalParameter(ECFParameter* parameter)
{
	parameters.push_back(parameter);
	return parameter;
}

void ECFObject::dropParameter(ECFParameter *parameter)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (parameters[i] == parameter) {
			parameters.erase(parameters.begin() + i);
			return;
		}
	}
	eslog::globalerror("ECF ERROR: Object cannot drop parameter '%s'\n", parameter->name.c_str());
}

void ECFObject::dropAllParameters()
{
	parameters.clear();
}

ECFParameter* ECFObject::addSeparator()
{
	registerParameter("SEPARATOR", new ECFSeparator(), ECFMetaData().setdatatype({ ECFDataType::SEPARATOR }));
	return parameters.back();
}

ECFParameter* ECFObject::addSpace()
{
	registerParameter("SPACE", new ECFSeparator(), ECFMetaData().setdatatype({ ECFDataType::SPACE }));
	return parameters.back();
}

void ECFObject::moveLastBefore(const std::string &name)
{
	ECFParameter *last = parameters.back();
	auto position = std::find(parameters.begin(), parameters.end(), getWithError(name));
	parameters.pop_back();
	parameters.insert(position, last);
}

ECFParameter* ECFObject::registerParameter(const std::string &name, ECFParameter *parameter, const ECFMetaData &metadata)
{
	registeredParameters.push_back(parameter);
	parameters.push_back(parameter);
	parameters.back()->name = name;
	parameters.back()->metadata = metadata;
	parameters.back()->defaultName();
	return parameters.back();
}

ECFParameter* ECFObject::registerPatternParameter(ECFParameter *parameter) const
{
	registeredParameters.push_back(parameter);
	parameter->name = metadata.pattern.front();
	parameter->metadata = metadata.suffix(1);
	parameter->defaultName();
	return parameter;
}

void ECFObject::forEachParameters(std::function<void(const ECFParameter*)> fnc, bool onlyAllowed, bool includePattern) const
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (onlyAllowed && !parameters[i]->metadata.isallowed()) {
			continue;
		}
		if (parameters[i]->isObject()) {
			fnc(parameters[i]);
			dynamic_cast<const ECFObject*>(parameters[i])->forEachParameters(fnc, onlyAllowed, includePattern);
		}
		if (parameters[i]->isValue()) {
			fnc(parameters[i]);
		}
	}
	
	ECFParameter* pattern;
	if (includePattern && (pattern = this->getPattern()) != NULL)
	{
		
		if (pattern->isObject()) {
			fnc(pattern);
			dynamic_cast<const ECFObject*>(pattern)->forEachParameters(fnc, onlyAllowed, includePattern);
		}
		if (pattern->isValue()) {
			fnc(pattern);
		}
	}
}

void ECFObject::forEachParameters(std::function<void(ECFParameter*)> fnc, bool onlyAllowed, bool includePattern, bool includeObjects)
{
	for (size_t i = 0; i < parameters.size(); i++) {
		if (onlyAllowed && !parameters[i]->metadata.isallowed()) {
			continue;
		}
		if (parameters[i]->isObject()) {
			if (includeObjects) { fnc(parameters[i]); }
			dynamic_cast<ECFObject*>(parameters[i])->forEachParameters(fnc, onlyAllowed, includePattern);
		}
		if (parameters[i]->isValue()) {
			fnc(parameters[i]);
		}
	}

	ECFParameter* pattern;
	if (includePattern && (pattern = this->getPattern()))
	{
		if (pattern->isObject()) {
			fnc(pattern);
			dynamic_cast<ECFObject*>(pattern)->forEachParameters(fnc, onlyAllowed, includePattern);
		}
		if (pattern->isValue()) {
			fnc(pattern);
		}
	}
}

void ECFObject::forEachParameters(
		std::function<void(const ECFParameter*, const std::string&)> fnc,
		const std::string& prefix, bool onlyAllowed, 
		bool includePattern
) const
{
	std::string new_prefix = prefix.empty() ? this->name : prefix + "." + this->name;
	for (size_t i = 0; i < parameters.size(); i++) {
		if (onlyAllowed && !parameters[i]->metadata.isallowed()) {
			continue;
		}
		fnc(parameters[i], new_prefix);
		if (parameters[i]->isObject()) {
			dynamic_cast<const ECFObject*>(parameters[i])->forEachParameters(fnc, new_prefix, onlyAllowed, includePattern);
		}
	}
	
	ECFParameter* pattern;
	if (includePattern && (pattern = this->getPattern()) != NULL)
	{
		if (pattern->isObject()) {
			pattern->name = this->name;
			fnc(pattern, prefix);
			dynamic_cast<const ECFObject*>(pattern)->forEachParameters(fnc, prefix, onlyAllowed, includePattern);
		}
		if (pattern->isValue()) {
			fnc(pattern, new_prefix);
		}
	}
}

ECFObject::~ECFObject()
{
	for (size_t i = 0; i < registeredParameters.size(); i++) {
		delete registeredParameters[i];
	}
}

