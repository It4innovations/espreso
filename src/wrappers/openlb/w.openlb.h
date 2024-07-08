
#ifndef SRC_WRAPPERS_OPENLB_W_OPENLB_H_
#define SRC_WRAPPERS_OPENLB_W_OPENLB_H_

namespace espreso {

struct OpenLBDataHolder;

class OpenLB {
public:
    static bool isLinked();

    OpenLB();
    ~OpenLB();

    bool set();
    bool update();
    bool solve();

protected:
    OpenLBDataHolder *external;
};

}

#endif /* SRC_WRAPPERS_OPENLB_W_OPENLB_H_ */
