
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "../loader.h"
#include "elements/elements.h"
#include "settings.h"

namespace espreso {
namespace input {

class Generator: public Loader {

protected:
	Generator(Mesh &mesh, const Settings &settings): Loader(mesh), _settings(settings) {};

	virtual ~Generator() { };

	void elements(std::vector<Element*> &elements)
	{
		elementsMesh(elements);
		elementsMaterials(elements);
	}

	virtual void elementsMesh(std::vector<Element*> &elements) = 0;
	virtual void elementsMaterials(std::vector<Element*> &elements) = 0;

	virtual void pickElementsInInterval(const std::vector<Element*> &elements, std::vector<Element*> &selection, const Interval &interval) =0;
	virtual void pickNodesInInterval(const std::vector<Element*> &nodes, std::vector<Element*> &selection, const Interval &interval) =0;

	virtual void generateFacesInInterval(std::vector<Element*> &faces, const Interval &interval) =0;
	virtual void generateEdgesInInterval(std::vector<Element*> &edges, const Interval &interval) =0;


	void materials(std::vector<Material> &materials);
	virtual void settings(
			std::vector<Evaluator*> &evaluators,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);

	void loadProperties(
			std::vector<Evaluator*> &evaluators,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes,
			const std::string &name,
			std::vector<std::string> parameters,
			std::vector<Property> properties);

	const Settings _settings;
};

}
}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
