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

#ifndef ASYNC_AS_SYNC_H
#define ASYNC_AS_SYNC_H

#include "Base.h"

namespace async
{

namespace as
{

/**
 * Asynchronous call via pthreads
 */
template<class Executor, typename InitParameter, typename Parameter>
class Sync : public Base<Executor, InitParameter, Parameter>
{
public:
	Sync()
	{
	}

	~Sync()
	{
	}

	unsigned int addSyncBuffer(const void* buffer, size_t size, bool clone = false)
	{
		return Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size, false);
	}

	unsigned int addBuffer(const void* buffer, size_t size, bool clone = false)
	{
		return Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size, buffer == 0L);
	}
	
	void resizeBuffer(unsigned int id, const void* buffer, size_t size)
	{
		Base<Executor, InitParameter, Parameter>::_resizeBuffer(id, buffer, size);
	}

	const void* buffer(unsigned int id) const
	{
		if (Base<Executor, InitParameter, Parameter>::origin(id))
			return Base<Executor, InitParameter, Parameter>::origin(id);

		return Base<Executor, InitParameter, Parameter>::_buffer(id);
	}

	void sendBuffer(unsigned int id, size_t size)
	{
	}

	/**
	 * Does nothing (call has already finished because it is synchronous)
	 */
	void wait()
	{
	}

	void call(const Parameter &parameters)
	{
		Base<Executor, InitParameter, Parameter>::call(parameters);
	}
};

}

}

#endif // ASYNC_AS_SYNC_H
