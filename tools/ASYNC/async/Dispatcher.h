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

#ifndef ASYNC_DISPATCHER_H
#define ASYNC_DISPATCHER_H

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#ifdef USE_ASYNC_MPI
#include "async/as/MPIScheduler.h"
#endif // USE_ASYNC_MPI

#include "Config.h"
#include "ModuleBase.h"

namespace async
{

class Dispatcher
{
private:
#ifdef USE_MPI
	async::as::MPIScheduler m_scheduler;

	MPI_Comm m_comm;
#endif // USE_MPI

	unsigned int m_groupSize;

public:
	Dispatcher() :
#ifdef USE_MPI
			m_comm(MPI_COMM_WORLD),
#endif // USE_MPI
			m_groupSize(Config::groupSize())
	{
	}

	~Dispatcher()
	{
		// Delete all modules so we can create a new dispatcher
		// probably only important for testing
		ModuleBase::modules().clear();
	}

#ifdef USE_MPI
	void setCommunicator(MPI_Comm comm)
	{
		m_comm = comm;
	}
#endif

	/**
	 * Use this to overwrite the group size set from the environment variable
	 *
	 * @param groupSize The group size (excl. the MPI executor)
	 * @deprecated Use {@link Config::setGroupSize}
	 */
	void setGroupSize(unsigned int groupSize)
	{
		if (Config::mode() == MPI)
			m_groupSize = groupSize;
	}

	/**
	 * Initialize the dispatcher
	 *
	 * This has to be called after {@link setCommunicator} and
	 * {@link setGroupSize}
	 */
	void init()
	{
#ifdef USE_MPI
		const std::vector<ModuleBase*>& modules = ModuleBase::modules();
		// Set the scheduler for all modules
		for (std::vector<ModuleBase*>::const_iterator i = modules.begin();
				i != modules.end(); ++i)
			(*i)->setScheduler(m_scheduler);

		if (Config::mode() == MPI)
			// Initialize the scheduler
			m_scheduler.setCommunicator(m_comm, m_groupSize);
#endif // USE_MPI
	}

	/**
	 * @return The groups size (or 1 for synchronous and asynchnchronous thread mode)
	 */
	unsigned int groupSize() const
	{
		return m_groupSize;
	}

#ifdef USE_MPI
	MPI_Comm groupComm() const
	{
		return m_scheduler.groupComm();
	}

	MPI_Comm commWorld() const
	{
		return m_scheduler.commWorld();
	}
#endif // USE_MPI

	/**
	 * @return True if the process is an MPI executor
	 */
	bool isExecutor() const
	{
#ifdef USE_MPI
		return m_scheduler.isExecutor();
#else // USE_MPI
		return false;
#endif // USE_MPI
	}

	/**
	 * This function will not return for MPI executors until all executors have been
	 * finalized. The function has to be called after all async {@link Module}s have
	 * been created.
	 *
	 * @return False if this rank is an MPI executor that does not contribute to the
	 *  computation.
	 */
	bool dispatch()
	{
#ifdef USE_MPI
		if (m_scheduler.isExecutor()) {
			const std::vector<ModuleBase*>& modules = ModuleBase::modules();
			// Initialize the executor modules
			for (std::vector<ModuleBase*>::const_iterator i = modules.begin();
					i != modules.end(); ++i)
				(*i)->setUp();

			// Run the executor loop
			m_scheduler.loop();

			// Finalize the executor modules
			for (std::vector<ModuleBase*>::const_iterator i = modules.begin();
					i != modules.end(); ++i)
				(*i)->tearDown();
			return false;
		}
#endif // USE_MPI

		return true;
	}

	void finalize()
	{
#ifdef USE_MPI
		m_scheduler.finalize();
#endif // USE_MPI
	}
};

}

#endif // ASYNC_DISPATCHER_H
