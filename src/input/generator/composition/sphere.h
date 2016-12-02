
#ifndef SRC_INPUT_GENERATOR_COMPOSITION_SPHERE_H_
#define SRC_INPUT_GENERATOR_COMPOSITION_SPHERE_H_

#include "../../loader.h"

#include "../primitives/triple.h"

namespace espreso {

struct SphereConfiguration;
enum class ELEMENT_TYPE;

namespace input {

struct BlockGenerator;

struct SphereSettings {

	SphereSettings();
	SphereSettings(const SphereConfiguration &configuration);

	ELEMENT_TYPE etype;

	double innerRadius, outerRadius;
	size_t clusters, layers;
	Triple<size_t> domains, elements;

	bool uniformDecomposition;
};

class Sphere: public Loader {

public:
	Sphere(Mesh &mesh, const SphereSettings &settings, size_t index, size_t size);
	virtual ~Sphere();

	static void load(const SphereConfiguration &configuration, Mesh &mesh, size_t index, size_t size);

	virtual void points(Coordinates &coordinates);
	virtual void elements(std::vector<Element*> &elements, std::vector<Element*> &faces, std::vector<Element*> &edges);
	virtual void materials(std::vector<Material> &materials);
	virtual void neighbours(std::vector<Element*> &nodes, std::vector<int> &neighbours);
	virtual void regions(
			std::vector<Evaluator*> &evaluators,
			std::vector<Region> &regions,
			std::vector<Element*> &elements,
			std::vector<Element*> &faces,
			std::vector<Element*> &edges,
			std::vector<Element*> &nodes);

	virtual bool partitiate(const std::vector<Element*> &nodes, std::vector<eslocal> &partsPtrs, std::vector<std::vector<Element*> > &fixPoints, std::vector<Element*> &corners);

protected:
	enum class SIDE : int {
		UP = 0,
		FRONT = 1,
		DOWN = 2,
		BACK = 3,
		LEFT = 4,
		RIGHT = 5
	};

	SphereSettings _settings;
	BlockGenerator* _block;
	Triple<size_t> _clusterOffset;
	Triple<size_t> _subnodes;

	size_t _index;
	size_t _size;

	SIDE _side;
	size_t _row, _col, _layer;
};


}
}


#endif /* SRC_INPUT_GENERATOR_COMPOSITION_SPHERE_H_ */
