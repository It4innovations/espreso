
#include "kernel.h"
#include "solverdataprovider/provider.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.h"
#include "config/holders/expression.h"
#include "mesh/store/nodestore.h"
#include "physics/system/builder/builder.h"
#include "wrappers/mpi/communication.h"

#include <cmath>
#include <algorithm>

using namespace espreso;

Kernel::InstanceFiller::InstanceFiller(esint nvectors)
: begin(0), end(0), interval(-1), invalid(0), DOFs(0), insertK(true), insertM(true), insertC(true), insertR(true), insertF(true),
  reduction(1), K(NULL), M(NULL), C(NULL), R(NULL), F(NULL), offset(NULL)
{
    Re.initVectors(nvectors);
    Fe.initVectors(nvectors);
}

void Kernel::smoothstep(double &smoothStep, double &derivation, double edge0, double edge1, double value, size_t order) const
{
    value = std::max(0.0, std::min((value - edge0) / (edge1 - edge0), 1.0));

    switch (order) {
    case 0:
        smoothStep = value;
        if (value == 0 || value == 1) {
            derivation = 0;
        } else {
            derivation = 1;
        }
        break;
    case 1:
        smoothStep = -2 * pow(value, 3) + 3 * pow(value, 2);
        derivation = -6 * pow(value, 2) + 6 * value;
        break;
    case 2:
        smoothStep = 6 * pow(value, 5) - 15 * pow(value, 4) + 10 * pow(value, 3);
        derivation = 30 * pow(value, 4) - 60*  pow(value, 3) + 30 * pow(value, 2);
        break;
    case 3:
        smoothStep = -20 * pow(value, 7) + 70 * pow(value, 6) - 84 * pow(value, 5) + 35 * pow(value, 4);
        derivation = -140 * pow(value, 6) + 420 * pow(value, 5) - 420 * pow(value, 4) + 140 * pow(value, 3);
        break;
    case 4:
        smoothStep = 70 * pow(value, 9) - 315 * pow(value, 8) + 540 * pow(value, 7) - 420 * pow(value, 6) + 126 * pow(value, 5);
        derivation = 630 * pow(value, 8) - 2520 * pow(value, 7) + 3780 * pow(value, 6) - 2520 * pow(value, 5) + 630 * pow(value, 4);
        break;
    case 5:
        smoothStep = -252 * pow(value, 11) + 1386 * pow(value, 10) - 3080 * pow(value, 9) + 3465 * pow(value, 8) - 1980 * pow(value, 7) + 462 * pow(value, 6);
        derivation = -2772 * pow(value, 10) + 13860 * pow(value, 9) - 27720 * pow(value, 8) + 27720 * pow(value, 7) - 13860 * pow(value, 6) + 2772 * pow(value, 5);
        break;
    case 6:
        smoothStep = 924 * pow(value, 13) - 6006 * pow(value, 12) + 16380 * pow(value, 11) - 24024 * pow(value, 10) + 20020 * pow(value, 9) - 9009 * pow(value, 8) + 1716 * pow(value, 7);
        derivation = 12012 * pow(value, 12) - 72072 * pow(value, 11) + 180180 * pow(value, 10) - 240240 * pow(value, 9) + 180180 * pow(value, 8) - 72072 * pow(value, 7) + 12012 * pow(value, 6);
        break;
    }

    derivation /= edge1 - edge0;
}


