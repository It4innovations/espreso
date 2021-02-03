
#include "distributed.composer.opt.h"

using namespace espreso;

void DistributedComposerOpt::assemble(const Builder &builder)
{
	if (kernel) {
		DistributedComposer::assemble(builder);
		return;
	}
	if (builder.matrices == Builder::Request::NONE) {
		return;
	}
}
