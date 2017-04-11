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

#ifndef ASYNC_EXECINFO_H
#define ASYNC_EXECINFO_H

#include <cassert>
#include <vector>

namespace async
{

/**
 * Buffer information send to the executor on each exec and execInit call
 */
class ExecInfo
{
private:
	/** The size for all buffers */
	std::vector<size_t> m_bufferSize;

public:
	virtual ~ExecInfo()
	{ }

	/**
	 * @return True, if this is an MPI executor
	 */
	virtual bool isExecutor() const
	{
		return false; // Default for sync and thread
	}

	unsigned int numBuffers() const
	{
		return m_bufferSize.size();
	}

	size_t bufferSize(unsigned int id) const
	{
		assert(id < numBuffers());
		return m_bufferSize[id];
	}

	/**
	 * @return Read-only pointer to the buffer (Useful for executors.)
	 */
	virtual const void* buffer(unsigned int id) const = 0;

protected:
	void _addBuffer(size_t size)
	{
		m_bufferSize.push_back(size);
	}

	void _removeBuffer(unsigned int id)
	{
		assert(id < numBuffers());
		m_bufferSize[id] = 0;
	}
};

}

#endif // ASYNC_EXECINFO_H