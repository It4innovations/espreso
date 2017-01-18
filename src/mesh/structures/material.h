
#ifndef SRC_MESH_MATERIALS_MATERIAL_H_
#define SRC_MESH_MATERIALS_MATERIAL_H_

#include <vector>
#include <string>
#include <algorithm>

#include "../../config/materialparameters.h"

namespace espreso {

class Coordinates;
class Evaluator;
struct Configuration;

class Material {

	friend std::ofstream& operator<<(std::ofstream& os, const Material &m);

public:
	Material(const Coordinates &coordinates, const Configuration &configuration);
	Material(std::ifstream &is, const Coordinates &coordinates): _coordinates(coordinates), _model(MATERIAL_MODEL::SIZE) {};

	virtual ~Material();

	const Evaluator* get(MATERIAL_PARAMETER parameter) const { return _values[static_cast<int>(parameter)]; }
	const Evaluator* get(size_t index) const { return _values[index]; }
	MATERIAL_MODEL getModel() const { return _model; }

	void set(size_t index, const std::string &value);
	void setModel(MATERIAL_MODEL model) { _model = model; }

protected:
	const Coordinates &_coordinates;

	MATERIAL_MODEL _model;
	std::vector<Evaluator*> _values;
};

}



#endif /* SRC_MESH_MATERIALS_MATERIAL_H_ */
