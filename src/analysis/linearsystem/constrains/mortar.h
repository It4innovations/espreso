
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
        int from , to;
        double value;
        Mortar(): pair(0), from(0), to(0), value(0) {}
        Mortar(int pair, int from, int to, double value): pair(pair), from(from), to(to), value(value) {}

        bool operator<(Mortar &other) { return pair == other.pair ? (from == other.from ? to < other.to : from < other.from) : pair < other.pair; }
        bool operator==(Mortar &other) { return pair == other.pair && from == other.from && to == other.to; }
        bool operator!=(Mortar &other) { return !(*this == other); }
    };

    struct Lambda {
        esint lambda, pair, id;
        Point normal;

        Lambda(): lambda(0), pair(0), id(0) {}
        Lambda(const Lambda &lambda) = default;
        Lambda(esint lambda, esint pair, esint id): lambda(lambda), pair(pair), id(id) {}
    };

    struct LambdaInfo: public Lambda {
        esint ndomains, doffset;

        LambdaInfo(const Lambda&lambda): Lambda(lambda), ndomains(0), doffset(0) {}
    };

    void assembleMortarInterface(std::vector<Mortar> &B);
    void synchronize(std::vector<Mortar> &B, std::vector<LambdaInfo> &lambdas, std::vector<int> &domains);

    esint ineq_begin, ineq_end;
    std::vector<LambdaInfo> lambdas;
    std::vector<int> domains;

    // pair, from, to, value -> normalized
    std::vector<Mortar> mortar;
};

}



#endif /* SRC_ANALYSIS_LINEARSYSTEM_CONSTRAINS_MORTAR_H_ */
