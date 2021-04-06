/**
 * @file
 *  This file is part of ASYNC
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2016-2017, Technische Universitaet Muenchen.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from this
 *     software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef ASYNC_AS_THREADBASE_H
#define ASYNC_AS_THREADBASE_H

#include <cstring>
#include <pthread.h>
#include <sched.h>

#include "utils/logger.h"

#include "Base.h"

namespace async
{

namespace as
{

/**
 * Base class for asynchronous calls via pthreads.
 *
 * This class is used by {@link Thread} and
 * {@link MPI}.
 */
template<class Executor, typename InitParameter, typename Parameter>
class ThreadBase : public Base<Executor, InitParameter, Parameter>
{
private:
	/**
	 * The possible phases, indicated whether we
	 * are betwenn wait() and call() or between
	 * call() and wait()
	 */
	enum Phase {
		/** We are between call() and wait() */
		EXEC_PHASE,
		/** We are between wait() and call() */
		SEND_PHASE
	};
	
	/** Async thread */
	pthread_t m_asyncThread;

	/** Mutex locked by the writer (caller) */
	pthread_spinlock_t m_writerLock;

	/** Mutex locked by the reader (callee) */
	pthread_mutex_t m_readerLock;

	/** Parameters for the next call */
	Parameter m_nextParams;

	/** The buffer we need to initialize */
	int m_initBuffer;
	
	/** The current phase */
	Phase m_phase;

	/** Mutex to wait for the buffer initialization to finish */
	pthread_spinlock_t m_initBufferLock;

	/** Shutdown the thread */
	bool m_shutdown;

protected:
	ThreadBase()
		: m_asyncThread(pthread_self()),
		  m_initBuffer(-1),
		  m_phase(EXEC_PHASE),
		  m_shutdown(false)
	{
		pthread_spin_init(&m_writerLock, PTHREAD_PROCESS_PRIVATE);
		pthread_mutex_init(&m_readerLock, 0L);
		pthread_spin_init(&m_initBufferLock, PTHREAD_PROCESS_PRIVATE);
	}

public:
	~ThreadBase()
	{
		finalize();

		pthread_spin_destroy(&m_initBufferLock);
		pthread_mutex_destroy(&m_readerLock);
		pthread_spin_destroy(&m_writerLock);
	}

	void setExecutor(Executor &executor)
	{
		Base<Executor, InitParameter, Parameter>::setExecutor(executor);

		// Lock the reader until data is available
		pthread_mutex_lock(&m_readerLock);

		// Lock writer until the memory is initialized
		pthread_spin_lock(&m_writerLock);

		// initBuffer is only unlocked when the buffer initialization is done
		pthread_spin_lock(&m_initBufferLock);

		if (pthread_create(&m_asyncThread, 0L, asyncThread, this) != 0)
			logError() << "ASYNC: Failed to start asynchronous thread";
	}
	
	void getAffinity(cpu_set_t &cpuSet)
	{
		pthread_getaffinity_np(m_asyncThread, sizeof(cpu_set_t), &cpuSet);
	}

	void setAffinity(const cpu_set_t &cpuSet)
	{
		pthread_setaffinity_np(m_asyncThread, sizeof(cpu_set_t), &cpuSet);
	}

	unsigned int addBuffer(const void* buffer, size_t size)
	{
		unsigned int id = Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size);

		// Now, initialize the buffer on the executor thread with zeros
		if (m_phase == EXEC_PHASE)
			pthread_spin_lock(&m_writerLock);
			
		m_initBuffer = id; // Mark for buffer fill
		pthread_mutex_unlock(&m_readerLock); // Similar to call() but without setting the parameters

		// Wait for the initialization to finish
		pthread_spin_lock(&m_initBufferLock);
		
		if (m_phase != EXEC_PHASE) // SEND_PHASE
			pthread_spin_lock(&m_writerLock);

		return id;
	}
	
	void resizeBuffer(unsigned int id, const void* buffer, size_t size)
	{
		Base<Executor, InitParameter, Parameter>::_resizeBuffer(id, buffer, size);
		
		assert(m_phase != EXEC_PHASE);
		
		// Initialize the buffer on the executor thread with zeros
		m_initBuffer = id; // Mark for buffer fill
		pthread_mutex_unlock(&m_readerLock); // Similar to call() but without setting the parameters

		// Wait for the initialization to finish
		pthread_spin_lock(&m_initBufferLock);
		
		pthread_spin_lock(&m_writerLock);
	}

	const void* buffer(unsigned int id) const
	{
		return Base<Executor, InitParameter, Parameter>::_buffer(id);
	}

	/**
	 * Wait for the asynchronous call to finish
	 */
	void wait()
	{
		pthread_spin_lock(&m_writerLock);
		m_phase = SEND_PHASE;
	}

	void call(const Parameter &parameters)
	{
		memcpy(&m_nextParams, &parameters, sizeof(Parameter));

		m_phase = EXEC_PHASE;
		pthread_mutex_unlock(&m_readerLock);
	}

	void finalize()
	{
		if (!Base<Executor, InitParameter, Parameter>::_finalize())
			return;

		// Shutdown the thread
		m_shutdown = true;
		pthread_mutex_unlock(&m_readerLock);
		pthread_join(m_asyncThread, 0L);
	}

private:
	/**
	 * Wrapper for the parent class because parent class function cannot be called directly
	 */
	void _call(const Parameter &parameters)
	{
		Base<Executor, InitParameter, Parameter>::call(parameters);
	}

private:
	static void* asyncThread(void* c)
	{
		ThreadBase* async = reinterpret_cast<ThreadBase*>(c);

		// Tell everyone that we are read to go
		pthread_spin_unlock(&async->m_writerLock);

		while (true) {
			// We assume that this lock happens before any unlock from the main thread
			pthread_mutex_lock(&async->m_readerLock);
			if (async->m_shutdown)
				break;

			if (async->m_initBuffer >= 0) {
				// Touch the memory on this thread
				unsigned int id = async->m_initBuffer;
				memset(async->_buffer(id), 0, async->bufferSize(id));
				async->m_initBuffer = -1;

				// Done
				pthread_spin_unlock(&async->m_initBufferLock);
			} else
				async->_call(async->m_nextParams);

			pthread_spin_unlock(&async->m_writerLock);
		}

		return 0L;
	}
};

}

}

#endif // ASYNC_AS_THREADBASE_H
