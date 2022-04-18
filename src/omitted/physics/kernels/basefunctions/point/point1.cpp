
#include "point1.h"
#include "basis/containers/serializededata.h"

using namespace espreso;

void Point1::setGaussPointsForOrder(int order)
{

}

void Point1::setBaseFunctions(Element &self)
{
	BaseFunctions::created(self);
}


