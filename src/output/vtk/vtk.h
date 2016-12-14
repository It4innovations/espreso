
#ifndef SRC_OUTPUT_VTK_VTK_H_
#define SRC_OUTPUT_VTK_VTK_H_

#include "../../assembler/constraints/constraints.h"
#include "../store.h"

namespace espreso {
namespace store {

class VTK: public Store {

public:
	int numb=1;
	VTK(const OutputConfiguration &output, const Mesh &mesh, const std::string &path);

	virtual void storeGeometry(size_t timeStep = -1);
	virtual void storeProperty(const std::string &name, const std::vector<Property> &properties, ElementType eType);
	virtual void storeValues(const std::string &name, size_t dimension, const std::vector<std::vector<double> > &values, ElementType eType);
	virtual void finalize();

	void store(std::vector<std::vector<double> > &displacement);

	static void gluing(const OutputConfiguration &output, const Mesh &mesh, const Constraints &constraints, const std::string &path, size_t dofs);

	static void mesh(const OutputConfiguration &output, const Mesh &mesh, const std::string &path, ElementType eType);
	static void fixPoints(const OutputConfiguration &output, const Mesh &mesh, const std::string &path);
	static void corners(const OutputConfiguration &output, const Mesh &mesh, const std::string &path);

protected:
	void computeCenters();
	Point shrink(const Point &p, size_t part) const;
	Point shrink(const Point &p, const Point &sCenter, const Point &cCenter) const;

	std::ofstream _os;
	std::vector<Point> _sCenters;
	Point _cCenter;
	ElementType _lastData;
};

}
}

#endif /* SRC_OUTPUT_VTK_VTK_H_ */
