
#include "w.pthread.h"
#include "basis/utilities/sysutils.h"

#include <pthread.h>

namespace espreso {

void* async(void *data);

struct ThreadControl {
	Pthread::Executor *executor;
	pthread_t thread;
	pthread_spinlock_t computationLock;
	pthread_mutex_t outputLock;

	bool finish;

	ThreadControl(Pthread::Executor *executor): executor(executor), finish(false)
	{
		pthread_spin_init(&computationLock, PTHREAD_PROCESS_PRIVATE);
		pthread_spin_lock(&computationLock);
		pthread_mutex_init(&outputLock, NULL);
		pthread_mutex_lock(&outputLock);

		pthread_attr_t attributes;
		pthread_attr_init(&attributes);

		cpu_set_t affinity;
		pthread_attr_getaffinity_np(&attributes, sizeof(cpu_set_t), &affinity);

		const int ncores = CPU_COUNT(&affinity);
		int core = ncores - 1;
		while (core >= 0) {
			if (CPU_ISSET(core, &affinity)) {
				core--;
			}
		}
		if (0 <= core) {
			cpu_set_t cpumask;
			CPU_ZERO(&cpumask);
			CPU_SET(core, &cpumask);

			pthread_attr_setaffinity_np(&attributes, sizeof(cpu_set_t), &cpumask);
		}

		pthread_create(&thread, &attributes, async, this);
		pthread_attr_destroy(&attributes);
	}

	~ThreadControl() {
		pthread_join(thread, NULL);
		pthread_spin_destroy(&computationLock);
		pthread_mutex_destroy(&outputLock);
	}
};

void* async(void *data)
{
	ThreadControl *threadControl = reinterpret_cast<ThreadControl*>(data);
	pthread_spin_unlock(&threadControl->computationLock);

	while (true) {
		pthread_mutex_lock(&threadControl->outputLock);
		if (threadControl->finish) {
			pthread_spin_unlock(&threadControl->computationLock);
			break;
		}
		threadControl->executor->call();
		pthread_spin_unlock(&threadControl->computationLock);
	}
	return NULL;
}

Pthread::Pthread(Executor *executor)
: _threadControl(new ThreadControl(executor))
{

}

Pthread::~Pthread()
{
	pthread_spin_lock(&_threadControl->computationLock);
	_threadControl->finish = true;
	pthread_mutex_unlock(&_threadControl->outputLock);
	pthread_spin_lock(&_threadControl->computationLock);
	delete _threadControl;
}

void Pthread::call()
{
	pthread_spin_lock(&_threadControl->computationLock);
	_threadControl->executor->copy();
	pthread_mutex_unlock(&_threadControl->outputLock);
}

}
