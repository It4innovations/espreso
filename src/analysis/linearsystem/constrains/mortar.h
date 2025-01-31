
#ifndef SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_MORTAR_H_
#define SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_MORTAR_H_

#include "esinfo/ecfinfo.h"
#include "feti/feti.h"

namespace espreso {

template <typename T>
struct MortarContact {
    void set(const step::Step &step, FETI<T> &feti);
    void update(const step::Step &step, FETI<T> &feti);

protected:
    struct Mortar {
        int pair;
        int from , to, domain;
        double value;
        Mortar(): pair(0), from(0), to(0), domain(-1), value(0) {}
        Mortar(int pair, int from, int to, double value): pair(pair), from(from), to(to), domain(-1), value(value) {}

        bool operator<(Mortar &other) {
            if (pair == other.pair) {
                if (from == other.from) {
                    if (to == other.to) {
                        return domain < other.domain;
                    }
                    return to < other.to;
                }
                return from < other.from;
            }
            return pair < other.pair;
        }
        bool operator==(Mortar &other) { return pair == other.pair && from == other.from && to == other.to && domain == other.domain; }
        bool operator!=(Mortar &other) { return !(*this == other); }
    };

    struct MortarInfo {
        int id, begin, end; // offset to Mortars
        Point normal;
        std::vector<int> domains;
    };

    void assembleMortarInterface(std::vector<Mortar> &B);
    void synchronize(FETI<T> &feti, std::vector<Mortar> &B);

    esint ineq_begin, ineq_end;

    // pair, from, to, value -> normalized
    std::vector<Mortar> mortar;
    std::map<int, MortarInfo> mInfo;
};

}



#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_MORTAR_H_ */
