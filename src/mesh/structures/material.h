
#ifndef MESH_STRUCTURES_MATERIAL_H_
#define MESH_STRUCTURES_MATERIAL_H_

#include "esbasis.h"
#include "../settings/evaluator.h"

namespace espreso {

class Material {

public:
	enum class MODEL {
			LINEAR_ELASTIC_ISOTROPIC = 0,
			LINEAR_ELASTIC_ORTHOTROPIC = 1,
			LINEAR_ELASTIC_ANISOTROPIC = 2
	};

	bool setParameter(const std::string &parameter, const std::string &value);

	MODEL model() const { return _model; }

	double density(size_t node) const { return _density->evaluate(node); }
	double termalCapacity(size_t node) const { return _termalCapacity->evaluate(node); }

	double youngModulusX(size_t node) const
	{
		return _youngModulus[0]->evaluate(node);
	}
	double youngModulusY(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _youngModulus[0]->evaluate(node);
		}
		return _youngModulus[1]->evaluate(node);
	}
	double youngModulusZ(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _youngModulus[0]->evaluate(node);
		}
		return _youngModulus[1]->evaluate(node);
	}

	double poissonRatioXY(size_t node) const
	{
			return _poissonRatio[0]->evaluate(node);
	}
	double poissonRatioXZ(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _poissonRatio[0]->evaluate(node);
		}
		return _poissonRatio[1]->evaluate(node);
	}
	double poissonRatioYZ(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _poissonRatio[0]->evaluate(node);
		}
		return _poissonRatio[2]->evaluate(node);
	}

	double termalExpansionX(size_t node) const
	{
		return _termalExpansion[0]->evaluate(node);
	}
	double termalExpansionY(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _termalExpansion[0]->evaluate(node);
		}
		return _termalExpansion[1]->evaluate(node);
	}
	double termalExpansionZ(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _termalExpansion[0]->evaluate(node);
		}
		return _termalExpansion[2]->evaluate(node);
	}

	double termalConductionX(size_t node) const
	{
		return _termalConduction[0]->evaluate(node);
	}
	double termalConductionY(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _termalConduction[0]->evaluate(node);
		}
		return _termalConduction[1]->evaluate(node);
	}
	double termalConductionZ(size_t node) const
	{
		if (_model == MODEL::LINEAR_ELASTIC_ISOTROPIC) {
			return _termalConduction[0]->evaluate(node);
		}
		return _termalConduction[2]->evaluate(node);
	}

	Material(const Coordinates &coordinates):
		_density(new ConstEvaluator(7850)), _termalCapacity(new ConstEvaluator(1)),
		_youngModulus{new ConstEvaluator(2.1e11), new ConstEvaluator(2.1e11), new ConstEvaluator(2.1e11)},
		_poissonRatio{new ConstEvaluator(0.3), new ConstEvaluator(0.3), new ConstEvaluator(0.3)},
		_termalExpansion{new ConstEvaluator(1), new ConstEvaluator(1), new ConstEvaluator(1)},
		_termalConduction{new ConstEvaluator(1), new ConstEvaluator(1), new ConstEvaluator(1)},
		_model(MODEL::LINEAR_ELASTIC_ISOTROPIC),
		_coordinates(&coordinates) {};

	Material(const Material &other)
	: _model(other._model), _coordinates(other._coordinates)
	{
		_density = other._density->copy();
		_termalCapacity = other._termalCapacity->copy();
		_youngModulus[0] = other._youngModulus[0]->copy();
		_youngModulus[1] = other._youngModulus[1]->copy();
		_youngModulus[2] = other._youngModulus[2]->copy();
		_poissonRatio[0] = other._poissonRatio[0]->copy();
		_poissonRatio[1] = other._poissonRatio[1]->copy();
		_poissonRatio[2] = other._poissonRatio[2]->copy();
		_termalExpansion[0] = other._termalExpansion[0]->copy();
		_termalExpansion[1] = other._termalExpansion[1]->copy();
		_termalExpansion[2] = other._termalExpansion[2]->copy();
		_termalConduction[0] = other._termalConduction[0]->copy();
		_termalConduction[1] = other._termalConduction[1]->copy();
		_termalConduction[2] = other._termalConduction[2]->copy();
	}

	Material& operator=(const Material &other)
	{
		if (this != &other) {
			delete _density;
			delete _termalCapacity;
			delete _youngModulus[0];
			delete _youngModulus[1];
			delete _youngModulus[2];
			delete _poissonRatio[0];
			delete _poissonRatio[1];
			delete _poissonRatio[2];
			delete _termalExpansion[0];
			delete _termalExpansion[1];
			delete _termalExpansion[2];
			delete _termalConduction[0];
			delete _termalConduction[1];
			delete _termalConduction[2];
			_model = other._model;
			_coordinates = other._coordinates;
			_density = other._density->copy();
			_termalCapacity = other._termalCapacity->copy();
			_youngModulus[0] = other._youngModulus[0]->copy();
			_youngModulus[1] = other._youngModulus[1]->copy();
			_youngModulus[2] = other._youngModulus[2]->copy();
			_poissonRatio[0] = other._poissonRatio[0]->copy();
			_poissonRatio[1] = other._poissonRatio[1]->copy();
			_poissonRatio[2] = other._poissonRatio[2]->copy();
			_termalExpansion[0] = other._termalExpansion[0]->copy();
			_termalExpansion[1] = other._termalExpansion[1]->copy();
			_termalExpansion[2] = other._termalExpansion[2]->copy();
			_termalConduction[0] = other._termalConduction[0]->copy();
			_termalConduction[1] = other._termalConduction[1]->copy();
			_termalConduction[2] = other._termalConduction[2]->copy();
		}
		return *this;
	}

	~Material()
	{
		delete _density;
		delete _termalCapacity;
		delete _youngModulus[0];
		delete _youngModulus[1];
		delete _youngModulus[2];
		delete _poissonRatio[0];
		delete _poissonRatio[1];
		delete _poissonRatio[2];
		delete _termalExpansion[0];
		delete _termalExpansion[1];
		delete _termalExpansion[2];
		delete _termalConduction[0];
		delete _termalConduction[1];
		delete _termalConduction[2];
	}

protected:
	Evaluator* _density;
	Evaluator* _termalCapacity;

	Evaluator* _youngModulus[3];
	Evaluator* _poissonRatio[3];
	Evaluator* _termalExpansion[3];
	Evaluator* _termalConduction[3];

	MODEL _model;
	const Coordinates* _coordinates;
};

}


#endif /* MESH_STRUCTURES_MATERIAL_H_ */
