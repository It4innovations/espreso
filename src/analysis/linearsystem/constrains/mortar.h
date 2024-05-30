
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_MORTAR_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_MORTAR_H_

#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct MortarContact {

    void set(const step::Step &step, FETI<T> &feti, int dofs);
    void update(const step::Step &step, FETI<T> &feti);

protected:
    struct Mortar {
        int pair;
        int from , to;
        double value;
        Mortar(): pair(0), from(0), to(0), value(0) {}
        Mortar(int pair, int from, int to, double value): pair(pair), from(from), to(to), value(value) {}

        bool operator<(Mortar &other) { return pair == other.pair ? (from == other.from ? to < other.to : from < other.from) : pair < other.pair; }
        bool operator==(Mortar &other) { return pair == other.pair && from == other.from && to == other.to; }
        bool operator!=(Mortar &other) { return !(*this == other); }
    };

    void assembleMortarInterface(std::vector<Mortar> &B);
    void synchronize(std::vector<Mortar> &B);

    std::vector<Mortar> mortar;
};

}



#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_MORTAR_H_ */
