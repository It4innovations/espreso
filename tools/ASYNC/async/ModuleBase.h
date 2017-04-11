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

#ifndef ASYNC_MODULEBASE_H
#define ASYNC_MODULEBASE_H

#include <vector>

#ifdef USE_MPI
#include "async/as/MPIScheduler.h"
#endif // USE_MPI

namespace async
{

class Dispatcher;

/**
 * Base class for asynchronous modules. Works closely together
 * with the {@link Dispatcher}.
 */
class ModuleBase
{
	friend class Dispatcher;
protected:
	ModuleBase()
	{
		modules().push_back(this);
	}

public:
	virtual ~ModuleBase()
	{ }

	/**
	 * Called at initialization. Is also called by the {@link Dispatcher}
	 * on MPI executors.
	 *
	 * Should at least call {@link setExecutor}(*this).
	 */
	virtual void setUp() = 0;

	/**
	 * Called after finalization. Is also called on MPI
	 * executors.
	 */
	virtual void tearDown()
	{ }

private:
	/**
	 * List of all I/O modules (required by the dispatcher)
	 */
	static std::vector<ModuleBase*>& modules()
	{
		// Use a function here to avoid an additional .cpp file
		static std::vector<ModuleBase*> moduleList;
		return moduleList;
	}

#ifdef USE_MPI
	/**
	 * Set the scheduler for this module.
	 */
	virtual void setScheduler(as::MPIScheduler &scheduler) = 0;
#endif // USE_MPI
};

}

#endif // ASYNC_MODULEBASE_H
