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

#include <cxxtest/TestSuite.h>

#include <cmath>

#include "async/as/Thread.h"
#include "Executor.h"

class TestThread : public CxxTest::TestSuite
{
private:
	pthread_spinlock_t m_lock;

	int m_value;

public:
	void setValue(int value)
	{
		m_value += value;
	}

	/**
	 * Tests the buffer alignement.
	 * Make sure that the environment variable ASYNC_BUFFER_ALIGNMENT=65536
	 * is set for this test.
	 */
	void testBuffer()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		int buffer = 42;
		async.addBuffer(&buffer, sizeof(int));

		async.wait();
		async.sendBuffer(0, sizeof(int));
		TS_ASSERT_EQUALS(*reinterpret_cast<const int*>(async.buffer(0)), 42);
		uintptr_t p = reinterpret_cast<uintptr_t>(async.buffer(0));
		TS_ASSERT_EQUALS(p % 65536, 0);
	}
};
