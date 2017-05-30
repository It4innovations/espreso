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

#ifndef ASYNC_AS_MPI_H
#define ASYNC_AS_MPI_H

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
class MPI : public MPIBase<Executor, InitParameter, Parameter>
{
public:
	MPI()
	{
	}

	~MPI()
	{
	}

	unsigned int addSyncBuffer(const void* buffer, size_t size, bool clone = false)
	{
		MPIBase<Executor, InitParameter, Parameter>::addBuffer(buffer, size, clone);
		return Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size, false);
	}

	unsigned int addBuffer(const void* buffer, size_t size, bool clone = false)
	{
		if (buffer == 0L)
			MPIBase<Executor, InitParameter, Parameter>::scheduler().addManagedBuffer(size);

		MPIBase<Executor, InitParameter, Parameter>::addBuffer(buffer, size, clone);
		return Base<Executor, InitParameter, Parameter>::_addBuffer(buffer, size, false);
	}
	
	void resizeBuffer(unsigned int id, const void* buffer, size_t size)
	{
		if (!Base<Executor, InitParameter, Parameter>::origin(id)) {
			// Resize the managed buffer
			MPIBase<Executor, InitParameter, Parameter>::scheduler().resizeManagedBuffer(
				Base<Executor, InitParameter, Parameter>::bufferSize(id), size);
		}
		
		MPIBase<Executor, InitParameter, Parameter>::resizeBuffer(id, buffer, size);
	}

	void removeBuffer(unsigned int id)
	{
		if (!Base<Executor, InitParameter, Parameter>::origin(id))
			MPIBase<Executor, InitParameter, Parameter>::scheduler().removeManagedBuffer(
				Base<Executor, InitParameter, Parameter>::bufferSize(id));

		MPIBase<Executor, InitParameter, Parameter>::removeBuffer(id);
	}

	void* managedBuffer(unsigned int id)
	{
		if (!Base<Executor, InitParameter, Parameter>::origin(id))
			return MPIBase<Executor, InitParameter, Parameter>::scheduler().managedBuffer();

		return 0L;
	}

	const void* buffer(unsigned int id) const
	{
		if (MPIBase<Executor, InitParameter, Parameter>::scheduler().isExecutor())
			return MPIBase<Executor, InitParameter, Parameter>::buffer(id);

		return Base<Executor, InitParameter, Parameter>::origin(id);
	}

	/**
	 * @warning Only the parameter from one task will be considered
	 */
	void call(const Parameter &parameters)
	{
		MPIBase<Executor, InitParameter, Parameter>::scheduler().sendParam(
			MPIBase<Executor, InitParameter, Parameter>::id(), parameters);
	}

private:
	bool useAsyncCopy() const
	{
		return false;
	}
};

}

}

#endif // ASYNC_AS_MPI_H
