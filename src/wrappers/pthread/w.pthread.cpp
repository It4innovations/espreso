
#include "w.pthread.h"
#include "basis/utilities/sysutils.h"

#include <pthread.h>
#include <atomic>

namespace espreso {

void* async(void *data);

struct SharedData {
	Pthread::Executor *executor;
	pthread_t thread;
	pthread_spinlock_t computationLock;
	pthread_mutex_t outputLock;

	std::atomic<int> tag;
	bool finish;

	SharedData(Pthread::Executor *executor): executor(executor), finish(false)
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
			pthread_spin_unlock(&shdata->computationLock);
			break;
		}
		int tag = std::atomic_load_explicit(&shdata->tag, std::memory_order_relaxed);
		std::atomic_thread_fence(std::memory_order_acquire);
		shdata->executor->call(tag);
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
	pthread_spin_lock(&_shdata->computationLock);
	delete _shdata;
}

void Pthread::execute(int tag)
{
	pthread_spin_lock(&_shdata->computationLock);
	_shdata->executor->copy(tag);
	std::atomic_thread_fence(std::memory_order_release);
	std::atomic_store_explicit(&_shdata->tag, tag, std::memory_order_relaxed);
	pthread_mutex_unlock(&_shdata->outputLock);
}

}
