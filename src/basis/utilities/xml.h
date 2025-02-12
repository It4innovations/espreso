
#ifndef SRC_BASIS_UTILITIES_XML_H_
#define SRC_BASIS_UTILITIES_XML_H_

#include <utility>
#include <string>
#include <vector>
#include <ostream>

namespace espreso {

class XML {
public:
    class Element {
    public:
        std::string name;
        std::string value;

        Element *parent;
        std::vector<std::pair<std::string, std::string> > attributes;
        std::vector<Element*> elements;

        XML::Element* attribute(const std::string &name, const std::string &value);
        XML::Element* element(const std::string &name);
        XML::Element* element(XML::Element* other);
        std::string ref(const std::string &name);

        std::pair<std::string, std::string>* getAttribute(const std::string &name);

        Element(): parent(NULL) {}
        Element(Element *parent): parent(parent) {}
    };

    std::string filename;
    std::string declaration;
    Element root;

    XML();
    ~XML();

    XML::Element* getFirstElement(const std::string &name, XML::Element *e = NULL);
    XML::Element* getNext(XML::Element *e);

    void print(std::ostream &os);
    void store(const std::string &filename);
    void load(const std::string &filename);
};

}



#endif /* SRC_BASIS_UTILITIES_XML_H_ */
