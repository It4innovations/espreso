
#ifndef SRC_WRAPPERS_PTHREAD_W_PTHREAD_H_
#define SRC_WRAPPERS_PTHREAD_W_PTHREAD_H_

namespace espreso {

struct ThreadControl;

class Pthread {
public:
    class Executor {
    public:
        virtual void call() = 0;
        virtual void copy() = 0;
        virtual ~Executor() {};
    };

    Pthread(Executor *executor);
    ~Pthread();

    void call();

protected:
    ThreadControl *_threadControl;
};

}



#endif /* SRC_WRAPPERS_PTHREAD_W_PTHREAD_H_ */
