
#ifndef SRC_MESH_STORE_CONTACTINFO_H_
#define SRC_MESH_STORE_CONTACTINFO_H_

namespace espreso {

struct ContactInfo {
    float gap, angle;
    bool self_contact;

    ContactInfo(): gap(0), angle(0), self_contact(false) {}
};

}

#endif /* SRC_MESH_STORE_CONTACTINFO_H_ */
