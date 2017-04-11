/**
 * @file
 *  This file is part of ASYNC
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2016, Technische Universitaet Muenchen.
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

#ifndef ASYNC_AS_MPIASYNC_H
#define ASYNC_AS_MPIASYNC_H

#include <mpi.h>

#include <cassert>
#include <vector>

#include "MPIBase.h"

namespace async
{

namespace as
{

/**
 * Asynchronous call via MPI
 */
template<class Executor, typename InitParameter, typename Parameter>
class MPIAsync : public MPIBase<Executor, InitParameter, Parameter>
{
private:
	/**
	 * Buffer description (on non-executors)
	 */
	struct BufInfo
	{
		/** A sychnronized buffer */
		bool sync;

		/**
		 * The number of asynchronous request required for this buffer.
		 *
		 * This is not counting the selecting isend.
		 */
		unsigned int requests;
	};

private:
	/** Buffer for the parameter */
	Parameter m_paramBuffer;

	/** Buffer description */
	std::vector<BufInfo> m_buffer;

	/** List of MPI requests */
	std::vector<MPI_Request> m_asyncRequests;

public:
	MPIAsync()
	{
		// One request always required for the parameters
		m_asyncRequests.push_back(MPI_REQUEST_NULL);
	}

	~MPIAsync()
	{
	}

	unsigned int addSyncBuffer(const void* buffer, size_t size, bool clone = false)
	{
		MPIBase<Executor, InitParameter, Parameter>::addBuffer(buffer, size, clone);
		unsigned int id = Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size, false);

		// We directly send sync buffers
		BufInfo bufInfo;
		bufInfo.sync = true;
		bufInfo.requests = 0;
		m_buffer.push_back(bufInfo);

		assert(m_buffer.size() == (Base<Executor, InitParameter, Parameter>::numBuffers()));

		return id;
	}

	/**
	 * @param bufferSize Should be 0 on the executor
	 */
	unsigned int addBuffer(const void* buffer, size_t size, bool clone = false)
	{
		MPIBase<Executor, InitParameter, Parameter>::addBuffer(buffer, size, clone, false);
		unsigned int id = Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size);

		// Initialize the requests
		unsigned int requests = 0;
		if (size > 0 && (!clone || MPIBase<Executor, InitParameter, Parameter>::scheduler().groupRank() == 0)) {
			requests = (size + MPIBase<Executor, InitParameter, Parameter>::maxSend() - 1)
				/ MPIBase<Executor, InitParameter, Parameter>::maxSend();
			m_asyncRequests.insert(m_asyncRequests.end(), requests*2, MPI_REQUEST_NULL);
		}

		BufInfo bufInfo;
		bufInfo.sync = false;
		bufInfo.requests = requests;
		m_buffer.push_back(bufInfo);

		assert(m_buffer.size() == (Base<Executor, InitParameter, Parameter>::numBuffers()));

		return id;
	}

	void removeBuffer(unsigned int id)
	{
		if (!m_buffer[id].sync) {
			m_asyncRequests.erase(m_asyncRequests.end()-m_buffer[id].requests*2,
				m_asyncRequests.end());
			m_buffer[id].requests = 0;
		}

		MPIBase<Executor, InitParameter, Parameter>::removeBuffer(id);
	}

	const void* buffer(unsigned int id) const
	{
		if (MPIBase<Executor, InitParameter, Parameter>::scheduler().isExecutor())
			return MPIBase<Executor, InitParameter, Parameter>::buffer(id);

		return Base<Executor, InitParameter, Parameter>::_buffer(id);
	}

	/**
	 * Wait for an asynchronous call to finish
	 */
	void wait()
	{
		// Wait for all requests first
		MPI_Waitall(m_asyncRequests.size(), &m_asyncRequests[0], MPI_STATUSES_IGNORE);

		// Wait for the call to finish
		MPIBase<Executor, InitParameter, Parameter>::wait();
	}

	/**
	 * @param id The id of the buffer
	 */
	void sendBuffer(unsigned int id, size_t size)
	{
		if (size == 0)
			return;

		assert(id < (Base<Executor, InitParameter, Parameter>::numBuffers()));

		if (m_buffer[id].sync) {
			MPIBase<Executor, InitParameter, Parameter>::sendBuffer(id, size);
			return;
		}

		// Only copy it to the local buffer
		assert((MPIBase<Executor, InitParameter, Parameter>::bufferPos(id)) + size
			<= (Base<Executor, InitParameter, Parameter>::bufferSize(id)));

		if (Base<Executor, InitParameter, Parameter>::origin(id)) {
			memcpy(Base<Executor, InitParameter, Parameter>::_buffer(id) +
					MPIBase<Executor, InitParameter, Parameter>::bufferPos(id),
				Base<Executor, InitParameter, Parameter>::origin(id) +
					MPIBase<Executor, InitParameter, Parameter>::bufferPos(id),
				size);
		}
		MPIBase<Executor, InitParameter, Parameter>::incBufferPos(id, size);
	}

	/**
	 * @warning Only the parameter from the last task will be considered
	 */
	void callInit(const InitParameter &parameters)
	{
		iSendAllBuffers();

		MPIBase<Executor, InitParameter, Parameter>::callInit(parameters);
	}

	/**
	 * @warning Only the parameter from the last task will be considered
	 */
	void call(const Parameter &parameters)
	{
		iSendAllBuffers();

		// Send parameters
		m_paramBuffer = parameters;
		MPIBase<Executor, InitParameter, Parameter>::scheduler().iSendParam(
			MPIBase<Executor, InitParameter, Parameter>::id(), m_paramBuffer);
	}

private:
	bool useAsyncCopy() const
	{
		return true;
	}

	/**
	 * Sends all buffers asynchronously
	 *
	 * Should only be used in asynchronous copy mode
	 */
	void iSendAllBuffers()
	{
		unsigned int nextRequest = 0;

		// Send all buffers
		for (unsigned int i = 0; i < Base<Executor, InitParameter, Parameter>::numBuffers(); i++) {
			size_t done = 0;
			for (unsigned int j = 0; j < m_buffer[i].requests; j++) {
				size_t send = std::min(MPIBase<Executor, InitParameter, Parameter>::maxSend(),
					MPIBase<Executor, InitParameter, Parameter>::bufferPos(i)-done);
				MPIRequest2 requests = MPIBase<Executor, InitParameter, Parameter>::scheduler().iSendBuffer(
					MPIBase<Executor, InitParameter, Parameter>::id(),
					i,
					Base<Executor, InitParameter, Parameter>::_buffer(i)+done,
					send);
				done += send;

				m_asyncRequests[nextRequest] = requests.r[0];
				m_asyncRequests[nextRequest+1] = requests.r[1];
				nextRequest += 2;
			}
		}

		assert(nextRequest == m_asyncRequests.size()-1);
	}
};

}

}

#endif // ASYNC_AS_MPIASYNC_H
