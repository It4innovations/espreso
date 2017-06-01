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

#ifndef ASYNC_MODULE_H
#define ASYNC_MODULE_H

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include "utils/env.h"

#ifdef USE_MPI
#include "async/as/MPI.h"
#include "async/as/MPIAsync.h"
#endif // USE_MPI
#include "async/as/Thread.h"
#include "async/as/Sync.h"

#include "Config.h"
#include "ModuleBase.h"

namespace async
{

template<class Executor, typename InitParameter, typename Parameter>
class Module : public ModuleBase
{
private:
	async::as::Base<Executor, InitParameter, Parameter>* m_async;

public:
	Module()
	{
		switch (Config::mode()) {
		case SYNC:
			m_async = new async::as::Sync<Executor, InitParameter, Parameter>();
			break;
		case THREAD:
			m_async = new async::as::Thread<Executor, InitParameter, Parameter>();
			break;
		case MPI:
#ifdef USE_MPI
			if (Config::useAsyncCopy())
				m_async = new async::as::MPIAsync<Executor, InitParameter, Parameter>();
			else
				m_async = new async::as::MPI<Executor, InitParameter, Parameter>();
#else // USE_MPI
			logError() << "Asynchronous MPI is not supported.";
#endif // USE_MPI
			break;
		}
	}

	virtual ~Module()
	{
		delete m_async;
	}

	void setExecutor(Executor &executor)
	{
		m_async->setExecutor(executor);
	}

	void init()
	{
		setUp();
	}

	unsigned int addSyncBuffer(const void* buffer, size_t size, bool clone = false)
	{
		return m_async->addSyncBuffer(buffer, size, clone);
	}

	unsigned int addBuffer(const void* buffer, size_t size, bool clone = false)
	{
		return m_async->addBuffer(buffer, size, clone);
	}
	
	void resizeBuffer(unsigned int id, const void* buffer, size_t size)
	{
		m_async->resizeBuffer(id, buffer, size);
	}

	void removeBuffer(unsigned int id)
	{
		m_async->removeBuffer(id);
	}

	unsigned int numBuffers() const
	{
		return m_async->numBuffers();
	}

	size_t bufferSize(unsigned int id) const
	{
		return m_async->bufferSize(id);
	}

	template<typename T>
	T managedBuffer(unsigned int id)
	{
		return static_cast<T>(m_async->managedBuffer(id));
	}

	const void* buffer(unsigned int id) const
	{
		return m_async->buffer(id);
	}

	/**
	 * Sends the complete buffer
	 */
	void sendBuffer(unsigned int id)
	{
		sendBuffer(id, bufferSize(id));
	}

	void sendBuffer(unsigned int id, size_t size)
	{
		m_async->sendBuffer(id, size);
	}

	void callInit(const InitParameter &param)
	{
		m_async->callInit(param);
	}

	void call(const Parameter &param)
	{
		m_async->call(param);
	}

	void wait()
	{
		m_async->wait();
	}

	void finalize()
	{
		m_async->finalize();

		tearDown();
	}

private:
#ifdef USE_MPI
	void setScheduler(as::MPIScheduler &scheduler)
	{
		m_async->setScheduler(scheduler);
	}
#endif // USE_MPI
};

}

#endif // ASYNC_MODULE_H
