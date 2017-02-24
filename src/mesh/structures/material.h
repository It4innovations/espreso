
#ifndef SRC_MESH_MATERIALS_MATERIAL_H_
#define SRC_MESH_MATERIALS_MATERIAL_H_

#include <vector>
#include <string>
#include <algorithm>

#include "../../configuration/material/parameters.h"
#include "../../configuration/physics.h"

namespace espreso {

class Coordinates;
class Evaluator;
struct Configuration;
struct CoordinateSystem;

struct MaterialCoordination {

	MaterialCoordination();
	MaterialCoordination(const Coordinates &coordinates, const CoordinateSystem &coordinateSystem);
	~MaterialCoordination();

	enum class Type {
		CARTESIAN,
		CYLINDRICAL,
		SPHERICAL
	};

	Type type;
	Evaluator* center[3];
	Evaluator* rotation[3];
};

class Material {

public:
	Material(const Coordinates &coordinates);
	Material(const Coordinates &coordinates, const Configuration &configuration, const CoordinateSystem &coordination);

	virtual ~Material();

	const Evaluator* get(MATERIAL_PARAMETER parameter) const { return _values[static_cast<int>(parameter)]; }
	MATERIAL_MODEL getModel(PHYSICS physics) const { return _models[(size_t)physics]; }

	void set(MATERIAL_PARAMETER parameter, const std::string &value);
	void set(MATERIAL_PARAMETER parameter, Evaluator* value);
	void setModel(PHYSICS physics, MATERIAL_MODEL model) { _models[(size_t)physics] = model; }

	const MaterialCoordination& coordination() const { return _coordination; }

	void store(std::ofstream& os);
	void load(std::ifstream& is);

protected:
	const Coordinates &_coordinates;

	std::vector<MATERIAL_MODEL> _models;
	std::vector<Evaluator*> _values;
	MaterialCoordination _coordination;
};

}



#endif /* SRC_MESH_MATERIALS_MATERIAL_H_ */
