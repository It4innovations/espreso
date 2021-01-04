
#include "w.pthread.h"
#include "basis/utilities/sysutils.h"

#include <pthread.h>

namespace espreso {

void* async(void *data);

struct SharedData {
	Pthread::Executor *executor;
	pthread_t thread;
	pthread_spinlock_t computationLock;
	pthread_mutex_t outputLock;

	int tag;
	bool finish;

	SharedData(Pthread::Executor *executor): executor(executor), finish(false)
	{
		pthread_spin_init(&computationLock, PTHREAD_PROCESS_PRIVATE);
		pthread_spin_lock(&computationLock);
		pthread_mutex_init(&outputLock, NULL);
		pthread_mutex_lock(&outputLock);

		cpu_set_t cpumask;
		CPU_ZERO(&cpumask);
		CPU_SET(utils::nprocs() - 1, &cpumask);

		pthread_attr_t attributes;
		pthread_attr_init(&attributes);
		pthread_attr_setaffinity_np(&attributes, sizeof(cpu_set_t), &cpumask);
		pthread_create(&thread, &attributes, async, this);
		pthread_attr_destroy(&attributes);
	}

	~SharedData() {
		pthread_join(thread, NULL);
		pthread_spin_destroy(&computationLock);
		pthread_mutex_destroy(&outputLock);
	}
};

void* async(void *data)
{
	SharedData *shdata = reinterpret_cast<SharedData*>(data);
	pthread_spin_unlock(&shdata->computationLock);

	while (true) {
		pthread_mutex_lock(&shdata->outputLock);
		if (shdata->finish) {
			break;
		}
		shdata->executor->call(shdata->tag);
		pthread_spin_unlock(&shdata->computationLock);
	}
	return NULL;
}

Pthread::Pthread(Executor *executor)
: _shdata(new SharedData(executor))
{

}

Pthread::~Pthread()
{
	pthread_spin_lock(&_shdata->computationLock);
	_shdata->finish = true;
	pthread_mutex_unlock(&_shdata->outputLock);
	delete _shdata;
}

void Pthread::call(int tag)
{
	pthread_spin_lock(&_shdata->computationLock);
	_shdata->tag = tag;
	pthread_mutex_unlock(&_shdata->outputLock);
}

}
