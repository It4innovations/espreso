
#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include "assembler.h"

class Dynamics: public Assembler {

public:
	Dynamics(const Instance &instance): Assembler(instance) { };

	void update();
	void solve();
};



#endif /* DYNAMICS_H_ */
