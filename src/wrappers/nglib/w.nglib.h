
#ifndef SRC_WRAPPERS_NGLIB_W_NGLIB_H_
#define SRC_WRAPPERS_NGLIB_W_NGLIB_H_

namespace espreso {

struct MeshData;

class NGLib {
public:
    static bool islinked();

    static void generate(MeshData &mesh);
};

}

#endif /* SRC_WRAPPERS_NGLIB_W_NGLIB_H_ */
