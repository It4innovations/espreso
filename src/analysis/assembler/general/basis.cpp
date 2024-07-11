
#include "basis.h"
#include "mesh/element.h"

namespace espreso {

double GaussPoints<Element::CODE::LINE2, 2, 2, 1>::w[];
double GaussPoints<Element::CODE::LINE2, 2, 2, 1>::N[];
double GaussPoints<Element::CODE::LINE2, 2, 2, 1>::dN[];
double GaussPoints<Element::CODE::LINE2, 2, 2, 1>::cw;
double GaussPoints<Element::CODE::LINE2, 2, 2, 1>::cN[];
double GaussPoints<Element::CODE::LINE2, 2, 2, 1>::cdN[];

double GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2>::w[];
double GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2>::N[];
double GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2>::dN[];
double GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2>::cw;
double GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2>::cN[];
double GaussPoints<Element::CODE::TRIANGLE3, 3, 6, 2>::cdN[];

double GaussPoints<Element::CODE::SQUARE4, 4, 4, 2>::w[];
double GaussPoints<Element::CODE::SQUARE4, 4, 4, 2>::N[];
double GaussPoints<Element::CODE::SQUARE4, 4, 4, 2>::dN[];
double GaussPoints<Element::CODE::SQUARE4, 4, 4, 2>::cw;
double GaussPoints<Element::CODE::SQUARE4, 4, 4, 2>::cN[];
double GaussPoints<Element::CODE::SQUARE4, 4, 4, 2>::cdN[];

double GaussPoints<Element::CODE::TETRA4, 4, 4, 3>::w[];
double GaussPoints<Element::CODE::TETRA4, 4, 4, 3>::N[];
double GaussPoints<Element::CODE::TETRA4, 4, 4, 3>::dN[];
double GaussPoints<Element::CODE::TETRA4, 4, 4, 3>::cw;
double GaussPoints<Element::CODE::TETRA4, 4, 4, 3>::cN[];
double GaussPoints<Element::CODE::TETRA4, 4, 4, 3>::cdN[];

double GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3>::w[];
double GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3>::N[];
double GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3>::dN[];
double GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3>::cw;
double GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3>::cN[];
double GaussPoints<Element::CODE::PYRAMID5, 5, 8, 3>::cdN[];

double GaussPoints<Element::CODE::PRISMA6, 6, 9, 3>::w[];
double GaussPoints<Element::CODE::PRISMA6, 6, 9, 3>::N[];
double GaussPoints<Element::CODE::PRISMA6, 6, 9, 3>::dN[];
double GaussPoints<Element::CODE::PRISMA6, 6, 9, 3>::cw;
double GaussPoints<Element::CODE::PRISMA6, 6, 9, 3>::cN[];
double GaussPoints<Element::CODE::PRISMA6, 6, 9, 3>::cdN[];

double GaussPoints<Element::CODE::HEXA8, 8, 8, 3>::w[];
double GaussPoints<Element::CODE::HEXA8, 8, 8, 3>::N[];
double GaussPoints<Element::CODE::HEXA8, 8, 8, 3>::dN[];
double GaussPoints<Element::CODE::HEXA8, 8, 8, 3>::cw;
double GaussPoints<Element::CODE::HEXA8, 8, 8, 3>::cN[];
double GaussPoints<Element::CODE::HEXA8, 8, 8, 3>::cdN[];

double GaussPoints<Element::CODE::LINE3, 3, 3, 1>::w[];
double GaussPoints<Element::CODE::LINE3, 3, 3, 1>::N[];
double GaussPoints<Element::CODE::LINE3, 3, 3, 1>::dN[];
double GaussPoints<Element::CODE::LINE3, 3, 3, 1>::cw;
double GaussPoints<Element::CODE::LINE3, 3, 3, 1>::cN[];
double GaussPoints<Element::CODE::LINE3, 3, 3, 1>::cdN[];

double GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2>::w[];
double GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2>::N[];
double GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2>::dN[];
double GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2>::cw;
double GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2>::cN[];
double GaussPoints<Element::CODE::TRIANGLE6, 6, 6, 2>::cdN[];

double GaussPoints<Element::CODE::SQUARE8, 8, 9, 2>::w[];
double GaussPoints<Element::CODE::SQUARE8, 8, 9, 2>::N[];
double GaussPoints<Element::CODE::SQUARE8, 8, 9, 2>::dN[];
double GaussPoints<Element::CODE::SQUARE8, 8, 9, 2>::cw;
double GaussPoints<Element::CODE::SQUARE8, 8, 9, 2>::cN[];
double GaussPoints<Element::CODE::SQUARE8, 8, 9, 2>::cdN[];

double GaussPoints<Element::CODE::TETRA10, 10, 15, 3>::w[];
double GaussPoints<Element::CODE::TETRA10, 10, 15, 3>::N[];
double GaussPoints<Element::CODE::TETRA10, 10, 15, 3>::dN[];
double GaussPoints<Element::CODE::TETRA10, 10, 15, 3>::cw;
double GaussPoints<Element::CODE::TETRA10, 10, 15, 3>::cN[];
double GaussPoints<Element::CODE::TETRA10, 10, 15, 3>::cdN[];

double GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3>::w[];
double GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3>::N[];
double GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3>::dN[];
double GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3>::cw;
double GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3>::cN[];
double GaussPoints<Element::CODE::PYRAMID13, 13, 14, 3>::cdN[];

double GaussPoints<Element::CODE::PRISMA15, 15, 9, 3>::w[];
double GaussPoints<Element::CODE::PRISMA15, 15, 9, 3>::N[];
double GaussPoints<Element::CODE::PRISMA15, 15, 9, 3>::dN[];
double GaussPoints<Element::CODE::PRISMA15, 15, 9, 3>::cw;
double GaussPoints<Element::CODE::PRISMA15, 15, 9, 3>::cN[];
double GaussPoints<Element::CODE::PRISMA15, 15, 9, 3>::cdN[];

double GaussPoints<Element::CODE::HEXA20, 20, 8, 3>::w[];
double GaussPoints<Element::CODE::HEXA20, 20, 8, 3>::N[];
double GaussPoints<Element::CODE::HEXA20, 20, 8, 3>::dN[];
double GaussPoints<Element::CODE::HEXA20, 20, 8, 3>::cw;
double GaussPoints<Element::CODE::HEXA20, 20, 8, 3>::cN[];
double GaussPoints<Element::CODE::HEXA20, 20, 8, 3>::cdN[];

}


