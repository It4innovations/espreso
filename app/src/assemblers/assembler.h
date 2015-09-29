
#ifndef ASSEMBLER_H_
#define ASSEMBLER_H_

#include "../instance.h"

class Assembler {

public:
	virtual ~Assembler() {};

	virtual void init() =0;
	virtual void pre_solve_update() =0;
	virtual void post_solve_update() =0;
	virtual void solve() =0;
	virtual void finalize() =0;

protected:
	Assembler(const Instance &instance): _instance(instance) { };

	const Instance &_instance;
};



#endif /* ASSEMBLER_H_ */
