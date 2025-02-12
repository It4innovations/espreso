
#include "prisma.h"

using namespace espreso;

bool Prisma::gpw(int order, std::vector<double> &r, std::vector<double> &s, std::vector<double> &t, std::vector<double> &w)
{
    switch (order) {
    case 0: break;
    default: return false;
    }
    return true;
}

int Prisma::maxorder()
{
    return 0;
}
