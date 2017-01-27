
#ifndef OUTPUT_ESDATA_ESDATA_H_
#define OUTPUT_ESDATA_ESDATA_H_

#include <cstdlib>
#include <string>

#include "../datastore.h"

namespace espreso {

class Mesh;
class Coordinates;

namespace store {

class ESPRESOBinaryFormat: public DataStore {

public:
	static void store(const Mesh &mesh, const std::string &path);

protected:
	ESPRESOBinaryFormat(const Mesh &mesh, const std::string &path);

	void coordinates();
	void elements();
	void materials();
	void regions();
	void boundaries();
};


}
}


#endif /* OUTPUT_ESDATA_ESDATA_H_ */
