
#ifndef SRC_CONFIG_READER_H_
#define SRC_CONFIG_READER_H_

#include "tokenizer.h"

namespace espreso {

class Reader {

public:
	static void read(const std::string &file) { _read(file, {}); }
	static void read(int* argc, char ***argv) { _read(argc, argv); }

	static void print();
	static void store();

private:
	static void _read(const std::string &file, const std::vector<std::string> &args);
	static void _read(int* argc, char ***argv);
};

}



#endif /* SRC_CONFIG_READER_H_ */
