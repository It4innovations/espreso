
#ifndef SRC_WRAPPERS_PAPI_W_PAPI_H_
#define SRC_WRAPPERS_PAPI_W_PAPI_H_

namespace espreso {

class PAPI {

public:
	PAPI();
	~PAPI();

	void init();
	long read();
	static bool isValid() { return status == 0; }

private:
	static int status;
	int set;
};

}

#endif /* SRC_WRAPPERS_PAPI_W_PAPI_H_ */
