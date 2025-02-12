
#ifndef SRC_OUTPUT_VISUALIZATION_WRITER_NETGENWRITTER_H_
#define SRC_OUTPUT_VISUALIZATION_WRITER_NETGENWRITTER_H_

#include "basis/io/outputfile.h"
#include "mesh/element.h"

namespace espreso {

struct NetgenASCIIWritter: public OutputFilePack {

    void int32(int value)
    {
        insert(snprintf(buffer, bsize, "%d", value));
    }

    void int32s(int value)
    {
        insert(snprintf(buffer, bsize, "%d ", value));
    }

    void int32ln(int value)
    {
        insert(snprintf(buffer, bsize, "%d\n", value));
    }

    void float32(float value)
    {
        insert(snprintf(buffer, bsize, "%f ", value));
    }

    void float32ln(float value)
    {
        insert(snprintf(buffer, bsize, "%f\n", value));
    }
};
}


#endif /* SRC_OUTPUT_VISUALIZATION_WRITER_NETGENWRITTER_H_ */
