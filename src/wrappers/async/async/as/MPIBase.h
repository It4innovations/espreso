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

#ifndef ASYNC_AS_MPIBASE_H
#define ASYNC_AS_MPIBASE_H

#include <mpi.h>

#include <cassert>
#include <cstring>
#include <vector>

#include "async/Config.h"
#include "ThreadBase.h"
#include "MPIScheduler.h"

namespace async
{

namespace as
{

/**
 * Asynchronous call via MPI
 *
 * @warning This class behaves very different depending on executor and non-executor
 *  ranks. Some variables are only available on non-executors while others are only
 *  available on executors.
 */
template<class Executor, typename InitParameter, typename Parameter>
class MPIBase : public ThreadBase<Executor, InitParameter, Parameter>, private Scheduled
{
private:
	struct BufInfo
	{

		/** True if this a clone buffer */
		bool clone;

		/** Current position of the buffer */
		size_t position;
	};

	/**
	 * Buffer description on the executor
	 */
	struct ExecutorBufInfo
	{
		/** True if the buffer uses synchronized sends */
		bool sync;

		/** Offsets for all tasks */
		unsigned long* offsets;

		/** Next writing position for all tasks */
		size_t* positions;

		/** Number of asychronous buffer chunks we receive for this buffer */
		unsigned int bufferChunks;
	};

private:
	/** The max amount that should be transfered in a single MPI send operation */
	const size_t m_maxSend;

	/** The scheduler */
	MPIScheduler* m_scheduler;

	/** The identifier for this module call */
	int m_id;

	/** Total number of asynchronous buffer chunks (only on the executor rank) */
	unsigned int m_numBufferChunks;

	/** Buffer description on non-executors */
	std::vector<BufInfo> m_buffer;

	/** Buffer description on executors */
	std::vector<ExecutorBufInfo> m_executorBuffer;

public:
	MPIBase()
		: m_maxSend(async::Config::maxSend()),
		  m_scheduler(0L),
		  m_id(-1),
		  m_numBufferChunks(0)
	{
	}

	~MPIBase()
	{
		finalize();
	}

	void setScheduler(MPIScheduler &scheduler)
	{
		m_scheduler = &scheduler;

		// Add this to the scheduler
		m_id = m_scheduler->addScheduled(this);
	}

	/**
	 * @param executor
	 */
	void setExecutor(Executor &executor)
	{
		// Initialization on the executor
		if (m_scheduler->isExecutor())
			ThreadBase<Executor, InitParameter, Parameter>::setExecutor(executor);
	}

	void removeBuffer(unsigned int id)
	{
		m_scheduler->removeBuffer(m_id, id);

		Base<Executor, InitParameter, Parameter>::removeBuffer(id);
	}

	const void* buffer(unsigned int id) const
	{
		if (m_scheduler->isExecutor())
			return ThreadBase<Executor, InitParameter, Parameter>::buffer(id);

		return 0L;
	}

	/**
	 * @param id The id of the buffer
	 */
	void sendBuffer(unsigned int id, size_t size)
	{
		if (size == 0)
			return;

		assert(id < (Base<Executor, InitParameter, Parameter>::numBuffers()));

		if (isClone(id) && m_scheduler->groupRank() != 0)
			return;

		const uint8_t* buffer = Base<Executor, InitParameter, Parameter>::origin(id);
		if (!buffer)
			buffer = m_scheduler->managedBuffer();

		// We need to send the buffer in 1 GB chunks
		for (size_t done = 0; done < size; done += maxSend()) {
			size_t send = std::min(maxSend(), size-done);

			m_scheduler->sendBuffer(m_id, id, buffer + bufferPos(id), send);
			incBufferPos(id, send);
		}
	}

	/**
	 * Wait for an asynchronous call to finish
	 */
	void wait()
	{
		// Wait for the call to finish
		m_scheduler->wait(m_id);

		resetBufferPosition();
	}

	/**
	 * @warning Only the parameter from one task will be considered
	 */
	void callInit(const InitParameter &parameters)
	{
		m_scheduler->sendInitParam(m_id, parameters);

		resetBufferPosition();
	}

	void finalize()
	{
		if (!Base<Executor, InitParameter, Parameter>::_finalize())
			return;

		if (m_id >= 0 && !m_scheduler->isExecutor())
			m_scheduler->sendFinalize(m_id);
	}

protected:
	size_t maxSend() const
	{
		return m_maxSend;
	}

	int id() const
	{
		return m_id;
	}

	MPIScheduler& scheduler()
	{
		return *m_scheduler;
	}

	const MPIScheduler& scheduler() const
	{
		return *m_scheduler;
	}

	unsigned int addBuffer(const void* buffer, size_t size, bool clone = false, bool sync = true)
	{
		assert(m_scheduler);
		assert(!m_scheduler->isExecutor());

		m_scheduler->addBuffer(m_id,
			Base<Executor, InitParameter, Parameter>::numBuffers());

		return _addBuffer(buffer, size, clone, sync);
	}
	
	void resizeBuffer(unsigned int id, const void* buffer, size_t size)
	{
		assert(m_scheduler);
		assert(!m_scheduler->isExecutor());
		
		// Resize the buffer on the compute node
		Base<Executor, InitParameter, Parameter>::_resizeBuffer(id, buffer, size);
		
		m_scheduler->resizeBuffer(m_id, id);
		
		_resizeBuffer(id, buffer, size);
	}

	bool isClone(unsigned int id) const
	{
		return m_buffer[id].clone;
	}

	/**
	 * The current position of a buffer on non-executors
	 */
	size_t bufferPos(unsigned int id) const
	{
		return m_buffer[id].position;
	}

	/**
	 * Increment the current buffer position on a non-executor
	 */
	void incBufferPos(unsigned int id, size_t increment)
	{
		m_buffer[id].position += increment;
	}

private:
	/**
	 * Reset the buffer position on non-executors
	 */
	void resetBufferPosition()
	{
		for (typename std::vector<BufInfo>::iterator it = m_buffer.begin();
				it != m_buffer.end(); it++)
			it->position = 0;
	}

	unsigned int paramSize() const
	{
		return std::max(sizeof(InitParameter), sizeof(Parameter));
	}

	virtual bool useAsyncCopy() const = 0;

	bool useAsyncCopy(unsigned int id) const
	{
		assert(id < m_executorBuffer.size());

		return useAsyncCopy() && !m_executorBuffer[id].sync;
	}

	unsigned int numBufferChunks() const
	{
		return m_numBufferChunks;
	}

	void _addBuffer(bool sync)
	{
		_addBuffer(0L, 0, false, sync);
	}

	unsigned int _addBuffer(const void* origin, unsigned long size, bool clone = false, bool sync = true)
	{
		// If this buffer is a clone, only the first rank will send it
		if (clone && m_scheduler->groupRank() != 0)
			size = 0;

		int executorRank = m_scheduler->groupSize()-1;

		// Compute buffer size and offsets
		unsigned long* bufferOffsets = 0L;
		if (m_scheduler->isExecutor()) {
			assert(origin == 0L);
			assert(size == 0);

			bufferOffsets = new unsigned long[m_scheduler->groupSize()];
		}

		MPI_Gather(&size, 1, MPI_UNSIGNED_LONG, bufferOffsets,
				1, MPI_UNSIGNED_LONG, executorRank, m_scheduler->privateGroupComm());

		if (m_scheduler->isExecutor()) {
			// Compute total size, offsets and buffer chunks
			size = 0;
			unsigned int bufferChunks = 0;
			for (int i = 0; i < m_scheduler->groupSize()-1; i++) {
				// Compute offsets from the size
				unsigned long bufSize = bufferOffsets[i];
				bufferOffsets[i] = size;

				if (bufSize > 0) {
					// Increment the total buffer size
					size += bufSize;

					if (!sync)
						// Increment the number of buffer chunks
						bufferChunks += (bufSize + m_maxSend - 1) / m_maxSend;
				}
			}

			m_numBufferChunks += bufferChunks;

			// Create the buffer
			ThreadBase<Executor, InitParameter, Parameter>::addBuffer(0L, size);

			ExecutorBufInfo executorBufInfo;
			executorBufInfo.sync = sync;
			executorBufInfo.offsets = bufferOffsets;
			executorBufInfo.bufferChunks = bufferChunks;

			// Initialize the current position
			executorBufInfo.positions = new size_t[m_scheduler->groupSize()-1];
			memset(executorBufInfo.positions, 0, (m_scheduler->groupSize()-1) * sizeof(size_t));

			m_executorBuffer.push_back(executorBufInfo);
		} else {
			BufInfo bufInfo;
			bufInfo.clone = clone;
			bufInfo.position = 0;

			m_buffer.push_back(bufInfo);
		}

		return Base<Executor, InitParameter, Parameter>::numBuffers()-1;
	}
	
	void _resizeBuffer(unsigned int id)
	{
		_resizeBuffer(id, 0L, 0);
	}
	
	void _resizeBuffer(unsigned int id, const void* buffer, size_t size)
	{
		if (!m_scheduler->isExecutor()) {
			if (m_buffer[id].clone && m_scheduler->groupRank() != 0)
				size = 0;
		}

		int executorRank = m_scheduler->groupSize()-1;

		// Compute buffer size and offsets
		unsigned long* bufferOffsets = 0L;
		if (m_scheduler->isExecutor()) {
			assert(size == 0);

			bufferOffsets = m_executorBuffer[id].offsets;
		}

		MPI_Gather(&size, 1, MPI_UNSIGNED_LONG, bufferOffsets,
				1, MPI_UNSIGNED_LONG, executorRank, m_scheduler->privateGroupComm());
		
		if (m_scheduler->isExecutor()) {
			// Compute total size, offsets and buffer chunks
			size = 0;
			unsigned int bufferChunks = 0;
			for (int i = 0; i < m_scheduler->groupSize()-1; i++) {
				// Compute offsets from the size
				unsigned long bufSize = bufferOffsets[i];
				bufferOffsets[i] = size;

				if (bufSize > 0) {
					// Increment the total buffer size
					size += bufSize;

					if (!m_executorBuffer[id].sync)
						// Increment the number of buffer chunks
						bufferChunks += (bufSize + m_maxSend - 1) / m_maxSend;
				}
			}

			m_numBufferChunks += bufferChunks - m_executorBuffer[id].bufferChunks;

			// Resize the buffer
			ThreadBase<Executor, InitParameter, Parameter>::resizeBuffer(id, 0L, size);

			m_executorBuffer[id].bufferChunks = bufferChunks;
		}
	}

	void _removeBuffer(unsigned int id)
	{
		assert(id < m_executorBuffer.size());

		m_numBufferChunks -= m_executorBuffer[id].bufferChunks;
		m_executorBuffer[id].bufferChunks = 0;

		delete [] m_executorBuffer[id].offsets;
		m_executorBuffer[id].offsets = 0L;
		delete [] m_executorBuffer[id].positions;
		m_executorBuffer[id].positions = 0L;

		Base<Executor, InitParameter, Parameter>::removeBuffer(id);
	}

	void* getBufferPos(unsigned int id, int rank, int size)
	{
		assert(rank < m_scheduler->groupSize()-1);
		assert(bufferOffset(id, rank)+size <= (Base<Executor, InitParameter, Parameter>::bufferSize(id)));

		void* buf = Base<Executor, InitParameter, Parameter>::_buffer(id)+bufferOffset(id, rank);
		m_executorBuffer[id].positions[rank] += size;
		return buf;
	}

	void _execInit(const void* paramBuffer)
	{
		const InitParameter* param = reinterpret_cast<const InitParameter*>(paramBuffer);
		Base<Executor, InitParameter, Parameter>::callInit(*param);

		resetBufferPositionOnEexecutor();
	}

	void _exec(const void* paramBuffer)
	{
		const Parameter* param = reinterpret_cast<const Parameter*>(paramBuffer);
		ThreadBase<Executor, InitParameter, Parameter>::call(*param);
	}

	void _wait()
	{
		ThreadBase<Executor, InitParameter, Parameter>::wait();

		resetBufferPositionOnEexecutor();
	}

	void _finalize()
	{
		for (typename std::vector<ExecutorBufInfo>::iterator it = m_executorBuffer.begin();
				it != m_executorBuffer.end(); ++it) {
			delete [] it->offsets;
			it->offsets = 0L;
			delete [] it->positions;
			it->positions = 0L;
		}

		ThreadBase<Executor, InitParameter, Parameter>::finalize();
	}

	/**
	 * @return The current offset on the buffer on the executor
	 */
	size_t bufferOffset(unsigned int id, int rank) const
	{
		return m_executorBuffer[id].offsets[rank]+m_executorBuffer[id].positions[rank];
	}

	/**
	 * Reset buffer positions
	 *
	 * Should only be called on the executor.
	 */
	void resetBufferPositionOnEexecutor()
	{
		for (typename std::vector<ExecutorBufInfo>::iterator it = m_executorBuffer.begin();
				it != m_executorBuffer.end(); it++) {
			if (it->positions)
				memset(it->positions, 0, (m_scheduler->groupSize()-1) * sizeof(size_t));
		}
	}
};

}

}

#endif // ASYNC_AS_MPIBASE_H
