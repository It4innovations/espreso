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

#include <mpi.h>

#include <cxxtest/TestSuite.h>

#include "async/as/MPIScheduler.h"

class TestMPIScheduler : public CxxTest::TestSuite
{
	int m_rank;

public:
	void setUp()
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
	}

	void testIsExecutor()
	{
		async::as::MPIScheduler scheduler;
		scheduler.setCommunicator(MPI_COMM_WORLD, 2);

		switch (m_rank) {
		case 2:
		case 4:
			TS_ASSERT(scheduler.isExecutor());
			break;
		default:
			TS_ASSERT(!scheduler.isExecutor());
		}
	}

	void testCommWorld()
	{
		async::as::MPIScheduler scheduler;
		scheduler.setCommunicator(MPI_COMM_WORLD, 2);

		int size;
		MPI_Comm_size(scheduler.commWorld(), &size);

		switch (m_rank) {
		case 2:
		case 4:
			TS_ASSERT_EQUALS(size, 2);
			break;
		default:
			TS_ASSERT_EQUALS(size, 3);
		}
	}

	void testGroupComm()
	{
		async::as::MPIScheduler scheduler;
		scheduler.setCommunicator(MPI_COMM_WORLD, 2);

		TS_ASSERT_EQUALS(scheduler.groupRank(), m_rank % 3);

		int size;

		switch (m_rank) {
		case 2:
		case 4:
			TS_ASSERT_EQUALS(scheduler.groupComm(), MPI_COMM_NULL);
			break;
		case 0:
		case 1:
			MPI_Comm_size(scheduler.groupComm(), &size);
			TS_ASSERT_EQUALS(size, 2);
			break;
		case 3:
			MPI_Comm_size(scheduler.groupComm(), &size);
			TS_ASSERT_EQUALS(size, 1);
		}

	}
};
