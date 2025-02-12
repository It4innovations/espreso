
#ifndef SRC_WRAPPERS_GMSH_W_GMSH_H_
#define SRC_WRAPPERS_GMSH_W_GMSH_H_

namespace espreso {

struct MeshData;

class GMSH {
public:
    static bool islinked();

    static void generate(MeshData &mesh);
};

}

#endif /* SRC_WRAPPERS_GMSH_W_GMSH_H_ */
