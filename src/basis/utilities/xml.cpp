
#include "xml.h"
#include "basis/utilities/parser.h"
#include "esinfo/eslog.hpp"
#include "esinfo/mpiinfo.h"

#include <fstream>
#include <functional>
#include <iostream>

using namespace espreso;

XML::XML()
: declaration("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
{

}

XML::~XML()
{
	std::function<void(Element*)> erase = [&] (Element *element) -> void {
		for (size_t i = 0; i < element->elements.size(); ++i) {
			erase(element->elements[i]);
			delete element->elements[i];
		}
	};

	erase(&root);
}

XML::Element* XML::getFirstElement(const std::string &name, XML::Element *e)
{
	std::function<XML::Element*(XML::Element*)> search = [&] (XML::Element *element) {
		if (StringCompare::caseInsensitiveEq(element->name, name)) {
			return element;
		}
		XML::Element *other = NULL;
		for (auto subelement = element->elements.begin(); subelement != element->elements.end(); ++subelement) {
			if ((other = search(*subelement)) != NULL) {
				return other;
			}
		}
		return other;
	};

	return search(e ? e : &root);
}

XML::Element* XML::getNext(XML::Element *e)
{
	bool after = false;
	std::function<XML::Element*(XML::Element*)> search = [&] (XML::Element *element) {
		if (StringCompare::caseInsensitiveEq(element->name, e->name)) {
			return element;
		}
		XML::Element *other = NULL;
		for (auto subelement = element->elements.begin(); subelement != element->elements.end(); ++subelement) {
			if ((other = search(*subelement)) != NULL) {
				if (after) {
					return other;
				}
				other = NULL;
				after = true;
			}
		}
		return other;
	};

	return search(&root);
}


XML::Element* XML::Element::attribute(const std::string &name, const std::string &value)
{
	attributes.push_back({ name, value});
	return this;
}

XML::Element* XML::Element::element(const std::string &name)
{
	elements.push_back(new Element(this));
	elements.back()->name = name;
	return elements.back();
}

XML::Element* XML::Element::element(XML::Element* other)
{
	elements.push_back(new Element(this));
	elements.back()->name = other->name;
	elements.back()->value = other->value;
	elements.back()->attributes = other->attributes;
	for (auto e = other->elements.begin(); e != other->elements.end(); ++e) {
		elements.back()->element(*e);
	}
	return elements.back();
}

std::string XML::Element::ref(const std::string &name)
{
	std::vector<XML::Element*> path = { this };
	while (path.back()->parent != NULL) {
		path.push_back(path.back()->parent);
	}

	std::string ret;
	for (auto it = path.rbegin(); it != path.rend(); ++it) {
		ret += "/" + (*it)->name;
	}
	ret += "[@Name=\"" + name + "\"]";

	return ret;
}

std::pair<std::string, std::string>* XML::Element::getAttribute(const std::string &name)
{
	for (auto it = attributes.begin(); it != attributes.end(); ++it) {
		if (StringCompare::caseInsensitiveEq(it->first, name)) {
			return &(*it);
		}
	}
	return NULL;
}

void XML::print(std::ostream &os)
{
	std::vector<Element*> stack = { &root };
	std::function<void(void)> _print = [&] () -> void {
		os << std::string(2 * stack.size() - 2, ' ') << "<" << stack.back()->name;
		for (auto attr = stack.back()->attributes.begin(); attr != stack.back()->attributes.end(); ++attr) {
			os << " " << attr->first << "=\"" << attr->second << "\"";
		}
		if (stack.back()->elements.size() == 0 && stack.back()->value.empty()) {
			os << "/>\n";
		} else {
			os << ">\n";
			if (!stack.back()->value.empty()) {
				os << std::string(2 * stack.size(), ' ') << stack.back()->value << "\n";
			}
			for (auto e = stack.back()->elements.begin(); e != stack.back()->elements.end(); ++e) {
				stack.push_back(*e);
				_print();
				stack.pop_back();
			}
			os << std::string(2 * stack.size() - 2, ' ') << "</" << stack.back()->name << ">\n";
		}
	};

	os << declaration << "\n";
	_print();
}

void XML::store(const std::string &filename)
{
	std::ofstream os(filename);
	print(os);
}

void XML::load(const std::string &filename)
{
	std::ifstream xml(filename);
	if (xml.fail() && info::mpi::rank == 0) {
		eslog::error("MESIO: XML file '%s' cannot be opened.\n", filename.c_str());
	}

	std::vector<char> buffer;
	XML::Element* e = NULL;

	auto isblank = [] (char c) {
		return c == ' ' || c == '\n' || c == '\t' || c == '\r';
	};

	char c;
	while (xml.get(c)) {
		if (c == '<') {
			if (buffer.size()) { // non-empty buffer is the value of the current element
				auto end = buffer.end() - 1;
				while (isblank(*end)) {
					--end;
				}
				e->value = std::string(buffer.begin(), end + 1);
				buffer.clear();
			}
			if (xml.peek() == '?') { // read declaration
				buffer.push_back('<');
				while (buffer.back() != '>') {
					buffer.push_back(xml.get());
				}
				declaration = std::string(buffer.begin(), buffer.end());
				buffer.clear();
				continue;
			}
			if (xml.peek() == '/') { // finish element tag
				while (c != '>') {
					xml.get(c);
				}
				if (e == NULL) {
					eslog::error("XML parser: incorrect XML (different number of open/close tags).\n");
				}
				e = e->parent;
				continue;
			}
			// new element
			if (e) {
				e = e->element("");
			} else {
				e = &root;
			}
			while (xml.get(c)) { // element name
				if (c == ' ' || c == '/' || c == '>') {
					e->name = std::string(buffer.begin(), buffer.end());
					buffer.clear();
					break;
				}
				buffer.push_back(c);
			}
			if (c == ' ') { // there are some parameters
				while (isblank(xml.peek())) {
					xml.get(c);
				}
				while (c != '>' && c != '/' && xml.peek() != '/') {
					std::string name;
					if (!isblank(c)) {
						buffer.push_back(c);
					}
					while (xml.get(c)) { // attribute name
						if (c == '=') {
							e->attribute(std::string(buffer.begin(), buffer.end()), "");
							buffer.clear();
							break;
						}
						if (!isblank(c)) {
							buffer.push_back(c);
						}
					}
					xml.get(c);
					xml.get(c);
					while (c != '"') {
						buffer.push_back(c);
						xml.get(c);
					}
					xml.get(c);
					e->attributes.back().second = std::string(buffer.begin(), buffer.end());
					buffer.clear();
					while (isblank(c)) {
						xml.get(c);
					}
				}
				while (c != '>') {
					if (c == '/') {
						e = e->parent;
					}
					xml.get(c);
				}
			}
		} else {
			if (buffer.size() || !isblank(c)) {
				buffer.push_back(c);
			}
		}
	}
}
