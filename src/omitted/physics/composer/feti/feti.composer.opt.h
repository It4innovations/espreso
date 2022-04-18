
#ifndef SRC_PHYSICS_COMPOSER_FETI_FETI_COMPOSER_OPT_H_
#define SRC_PHYSICS_COMPOSER_FETI_FETI_COMPOSER_OPT_H_

#include "feti.composer.h"

namespace espreso {

class FETIComposerOpt: public FETIComposer {

public:
	using FETIComposer::FETIComposer;

	void assemble(const Builder &builder);
};

}


#endif /* SRC_PHYSICS_COMPOSER_FETI_FETI_COMPOSER_OPT_H_ */
