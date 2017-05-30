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

#ifndef ASYNC_AS_THREAD_H
#define ASYNC_AS_THREAD_H

#include <cassert>
#include <cstring>
#include <sched.h>
#include <sys/sysinfo.h>

#include "async/Config.h"
#include "ThreadBase.h"

namespace async
{

namespace as
{

/**
 * Asynchronous call via pthreads
 */
template<class Executor, typename InitParameter, typename Parameter>
class Thread : public ThreadBase<Executor, InitParameter, Parameter>
{
private:
	/**
	 * Buffer description
	 */
	struct BufInfo
	{
		/** True if this is an init buffer */
		bool init;

		/** Current position of the buffer */
		size_t position;
	};

private:
	/** Buffer information */
	std::vector<BufInfo> m_buffer;

public:
	Thread()
	{
	}

	~Thread()
	{
	}

	void setExecutor(Executor &executor)
	{
		ThreadBase<Executor, InitParameter, Parameter>::setExecutor(executor);

		cpu_set_t oldCpuMask;
		ThreadBase<Executor, InitParameter, Parameter>::getAffinity(oldCpuMask);
		const int numCores = CPU_COUNT(&oldCpuMask); // Number of cores we have available

		int core = async::Config::getPinCore();
		if (core < 0)
			core = numCores + core;

		if (core < 0 || core >= numCores) {
			logWarning() << "Skipping async thread pining, invalid core" << core << "specified. Available cores:" << numCores;
			return;
		}
		
		// Get the real core
		const int totalCores = get_nprocs();
		int realCore = -1;
		while (core >= 0) {
			realCore++;
			if (CPU_ISSET(realCore, &oldCpuMask))
				core--;
			
			if (realCore >= totalCores)
				logError() << "Pinning failed. Not enough cores available.";
		}
		
		logDebug() << "Pinning executor to core" << realCore;
		
		cpu_set_t cpuMask;
		CPU_ZERO(&cpuMask);
		CPU_SET(realCore, &cpuMask);
		
		ThreadBase<Executor, InitParameter, Parameter>::setAffinity(cpuMask);
	}

	unsigned int addSyncBuffer(const void* buffer, size_t size, bool clone = false)
	{
		unsigned int id = Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size, false);
		BufInfo bufInfo;
		bufInfo.init = true;
		bufInfo.position = 0;
		m_buffer.push_back(bufInfo);

		assert(m_buffer.size() == (Base<Executor, InitParameter, Parameter>::numBuffers()));

		return id;
	}

	unsigned int addBuffer(const void* buffer, size_t size, bool clone = false)
	{
		unsigned int id = ThreadBase<Executor, InitParameter, Parameter>::addBuffer(buffer, size);
		BufInfo bufInfo;
		bufInfo.init = false;
		bufInfo.position = 0;
		m_buffer.push_back(bufInfo);

		assert(m_buffer.size() == (Base<Executor, InitParameter, Parameter>::numBuffers()));

		return id;
	}

	const void* buffer(unsigned int id) const
	{
		assert(id < m_buffer.size());

		if (m_buffer[id].init)
			return Base<Executor, InitParameter, Parameter>::origin(id);

		return ThreadBase<Executor, InitParameter, Parameter>::buffer(id);
	}

	void sendBuffer(unsigned int id, size_t size)
	{
		assert(id < (Base<Executor, InitParameter, Parameter>::numBuffers()));

		if (m_buffer[id].init)
			return;

		assert((Base<Executor, InitParameter, Parameter>::_buffer(id)));

		if (Base<Executor, InitParameter, Parameter>::origin(id)) {
			assert(m_buffer[id].position+size <= (Base<Executor, InitParameter, Parameter>::bufferSize(id)));

			memcpy(Base<Executor, InitParameter, Parameter>::_buffer(id)+m_buffer[id].position,
				Base<Executor, InitParameter, Parameter>::origin(id)+m_buffer[id].position,
				size);
		}
		m_buffer[id].position += size;
	}

	void wait()
	{
		ThreadBase<Executor, InitParameter, Parameter>::wait();

		resetBufferPosition();
	}

	void callInit(const InitParameter &parameters)
	{
		Base<Executor, InitParameter, Parameter>::callInit(parameters);

		resetBufferPosition();
	}

	void call(const Parameter &parameters)
	{
		ThreadBase<Executor, InitParameter, Parameter>::call(parameters);
	}

private:
	void resetBufferPosition()
	{
		for (typename std::vector<BufInfo>::iterator it = m_buffer.begin();
			it != m_buffer.end(); ++it)
			it->position = 0;
	}
};

}

}

#endif // ASYNC_AS_THREAD_H
