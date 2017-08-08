
#ifndef SRC_OUTPUT_DATASTORE_ESPRESOBINARYFORMAT_H_
#define SRC_OUTPUT_DATASTORE_ESPRESOBINARYFORMAT_H_

#include <cstdlib>
#include <string>

namespace espreso {

class Mesh;

class ESPRESOBinaryFormat {

public:
	static void prepareDirectories(const std::string &path, size_t parts);
	static void store(const Mesh &mesh, const std::string &path);

protected:
	ESPRESOBinaryFormat(const Mesh &mesh, const std::string &path);

	void metafile();
	void coordinates();
	void elements();
	void materials();
	void regions();
	void boundaries();

	const Mesh &_mesh;
	std::string _path;
};


}


#endif /* SRC_OUTPUT_DATASTORE_ESPRESOBINARYFORMAT_H_ */
