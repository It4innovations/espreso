
#ifndef SRC_WRAPPERS_PTHREAD_W_PTHREAD_H_
#define SRC_WRAPPERS_PTHREAD_W_PTHREAD_H_

namespace espreso {

struct SharedData;

class Pthread {
public:
	class Executor {
	public:
		virtual void call(int tag) = 0;
		virtual void copy(int tag) = 0;
		virtual ~Executor() {};
	};

	Pthread(Executor *executor);
	~Pthread();

	void execute(int tag);

protected:
	SharedData * _shdata;
};

class Async: public Pthread, public Pthread::Executor {
public:
	Async(Pthread::Executor *executor): Pthread(executor) {}
};

}



#endif /* SRC_WRAPPERS_PTHREAD_W_PTHREAD_H_ */
