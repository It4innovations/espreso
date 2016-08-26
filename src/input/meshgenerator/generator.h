
#ifndef INPUT_MESHGENERATOR_GENERATOR_H_
#define INPUT_MESHGENERATOR_GENERATOR_H_

#include "../loader.h"
#include "elements/elements.h"
#include "settings.h"

namespace espreso {
namespace input {

class Generator: public Loader {

protected:
	Generator(Mesh &mesh, const Settings &settings): Loader(mesh), _settings(settings)
	{
		switch (_settings.assembler) {
		case LinearElasticity:
			config::assembler::physics = config::assembler::LinearElasticity;
			config::assembler::timeSteps = 1;
			_DOFs = 3;
			break;
		case Temperature:
			config::assembler::physics = config::assembler::Temperature;
			config::assembler::timeSteps = 1;
			_DOFs = 1;
			break;
		case TransientElasticity:
			config::assembler::physics = config::assembler::TransientElasticity;
			_DOFs = 3;
			break;
		default:
			ESINFO(ERROR) << "Unknown assembler: ASSEMBLER = " << _settings.assembler;
		}
	};

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
	/*{
		auto set = [] (double &value, const std::map<std::string, double> &settings, const std::string &param) {
			if (settings.find(param) != settings.end()) {
				value = settings.find(param)->second;
			}
		};

		auto fillMaterial = [ &set ] (Material &material, const std::map<std::string, double> &settings) {
			set(material.model             , settings, "");
			set(material.density           , settings, "DENS");
			set(material.poissonRatio.z    , settings, "MIXY");
			set(material.poissonRatio.y    , settings, "MIXZ");
			set(material.poissonRatio.x    , settings, "MIYZ");
			set(material.youngModulus.x    , settings, "EX");
			set(material.youngModulus.y    , settings, "EY");
			set(material.youngModulus.z    , settings, "EZ");
			set(material.termalCapacity    , settings, "CP");
			set(material.termalConduction.x, settings, "KX");
			set(material.termalConduction.y, settings, "KY");
			set(material.termalConduction.z, settings, "KZ");
			set(material.termalExpansion.x , settings, "ALPHAX");
			set(material.termalExpansion.y , settings, "ALPHAY");
			set(material.termalExpansion.z , settings, "ALPHAZ");
		};

		materials.resize(2);
		fillMaterial(materials[0], _settings.material1);
		fillMaterial(materials[1], _settings.material2);
	}*/

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
	size_t _DOFs;
};

}
}


#endif /* INPUT_MESHGENERATOR_GENERATOR_H_ */
