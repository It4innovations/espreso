
#include "mpi.h"

#include "region.h"
#include "../elements/element.h"
#include "../structures/coordinates.h"
#include "../../basis/matrices/denseMatrix.h"

using namespace espreso;

Region::Region(ElementType eType)
: eType(eType), area(0), _destroy(true)
{
	_elements = new std::vector<Element*>();
}

Region::Region(ElementType eType, std::vector<Element*> &element)
: eType(eType), area(0), _elements(&element), _destroy(false)
{

}

void Region::computeArea(const Coordinates &coordinates) const
{
	double A = 0;
	for (size_t e = 0; e < elements().size(); e++) {

		DenseMatrix coords(elements()[e]->nodes(), 3), dND(1, 3);

		const std::vector<DenseMatrix> &dN = elements()[e]->dN();
		const std::vector<double> &weighFactor = elements()[e]->weighFactor();

		for (size_t n = 0; n < elements()[e]->nodes(); n++) {
			coords(n, 0) = coordinates[elements()[e]->node(n)].x;
			coords(n, 1) = coordinates[elements()[e]->node(n)].y;
			coords(n, 2) = coordinates[elements()[e]->node(n)].z;
		}

		if (elements()[e]->type() == Element::Type::LINE) {
			for (size_t gp = 0; gp < elements()[e]->gaussePoints(); gp++) {
				dND.multiply(dN[gp], coords);
				A += dND.norm() * weighFactor[gp];
			}
		}
		if (elements()[e]->type() == Element::Type::PLANE) {
			for (size_t gp = 0; gp < elements()[e]->gaussePoints(); gp++) {
				dND.multiply(dN[gp], coords);
				Point v2(dND(0, 0), dND(0, 1), dND(0, 2));
				Point v1(dND(1, 0), dND(1, 1), dND(1, 2));
				Point va = Point::cross(v1, v2);
				A += va.norm() * weighFactor[gp];
			}
		}
	}

	MPI_Allreduce(&A, &area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

Region::~Region()
{
	if (_destroy) {
		delete _elements;
	}
}

