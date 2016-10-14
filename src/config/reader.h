
#ifndef SRC_CONFIG_READER_H_
#define SRC_CONFIG_READER_H_

#include "tokenizer.h"

namespace espreso {

class Reader {

public:
	void read(int* argc, char ***argv) { _read(Tokenizer("espreso.config.new")); }

private:
	void _read(Tokenizer &&tokenizer);
};

}



#endif /* SRC_CONFIG_READER_H_ */
