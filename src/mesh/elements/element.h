#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include <cstring>

#include "esbasis.h"
#include "../settings/setting.h"

namespace espreso {

class Mesh;

class Element
{
	friend class Mesh;

public:
	enum Params {
		MATERIAL,
		CONSTANT,
		COORDINATES,
		BODY,
		PARAMS_SIZE
	};

	inline static bool match(const eslocal *indices, eslocal x, eslocal y)
	{
		return indices[x] == indices[y];
	}

	friend std::ofstream& operator<<(std::ofstream& os, const Element &e);

	inline bool operator==(const Element& other)
	{
		if (nodes() != other.nodes()) {
			return false;
		}
		return std::is_permutation(indices(), indices() + nodes(), other.indices());
	}

	virtual ~Element() {};

	virtual const std::vector<DenseMatrix>& dN() const = 0;
	virtual const std::vector<DenseMatrix>& N() const = 0;
	virtual const std::vector<double>& weighFactor() const = 0;

	virtual const std::vector<Property>& elementDOFs() const = 0;
	virtual const std::vector<Property>& faceDOFs() const = 0;
	virtual const std::vector<Property>& edgeDOFs() const = 0;
	virtual const std::vector<Property>& pointDOFs() const = 0;
	virtual const std::vector<Property>& midPointDOFs() const = 0;

	virtual eslocal nCommon() const = 0;
	virtual eslocal vtkCode() const = 0;

	virtual size_t faces() const = 0;
	virtual size_t edges() const = 0;
	virtual size_t nodes() const = 0;
	virtual size_t coarseNodes() const = 0;
	virtual size_t gaussePoints() const = 0;

	virtual Element* face(size_t index) const = 0;
	virtual Element* edge(size_t index) const = 0;
	eslocal& node(size_t index) { return indices()[index]; }
	const eslocal& node(size_t index) const { return indices()[index]; }

	virtual eslocal param(Params param) const =0;
	virtual void param(Params param, eslocal value) =0;

	Settings& settings() { return _settings; }
	const Settings& settings() const { return _settings; }

	std::vector<Evaluator*>& settings(Property property) { return _settings[property]; }
	const std::vector<Evaluator*>& settings(Property property) const { return _settings[property]; }

	std::vector<eslocal>& domains() { return _domains; }
	const std::vector<eslocal>& domains() const { return _domains; }

	std::vector<eslocal>& clusters() { return _clusters; }
	const std::vector<eslocal>& clusters() const { return _clusters; }

protected:
	virtual eslocal* indices() = 0;
	virtual const eslocal* indices() const = 0;
	virtual std::vector<eslocal> getNeighbours(size_t nodeIndex) const = 0;

	Settings _settings;
	std::vector<eslocal> _domains;
	std::vector<eslocal> _clusters;
};

inline std::ofstream& espreso::operator<<(std::ofstream& os, const Element &e)
{
	eslocal value = e.vtkCode();
	os.write(reinterpret_cast<const char *>(&value), sizeof(eslocal));
	os.write(reinterpret_cast<const char *>(e.indices()), sizeof(eslocal) * e.nodes());
	return os;
}

}


#endif /* ELEMENT_H_ */
