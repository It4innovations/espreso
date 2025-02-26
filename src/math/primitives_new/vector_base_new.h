
#ifndef SRC_MATH_PRIMITIVES_NEW_VECTOR_BASE_NEW_H_
#define SRC_MATH_PRIMITIVES_NEW_VECTOR_BASE_NEW_H_

#include <cstddef>



namespace espreso {



class VectorBase_new
{
public: // the user promises not to modify these values (I don't want to implement getters everywhere)
    size_t size = 0;
public:
    VectorBase_new() = default;
    VectorBase_new(const VectorBase_new &) = default;
    VectorBase_new(VectorBase_new &&) = default;
    VectorBase_new & operator=(const VectorBase_new &) = default;
    VectorBase_new & operator=(VectorBase_new &&) = default;
    virtual ~VectorBase_new() = default;
};



}



#endif /* SRC_MATH_PRIMITIVES_NEW_VECTOR_BASE_NEW_H_ */
