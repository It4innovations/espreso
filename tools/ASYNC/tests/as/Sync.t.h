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

#include "async/as/Sync.h"
#include "Executor.h"

class Param;
class TestSync : public CxxTest::TestSuite
{
private:
	int m_value;

public:
	void setValue(int value)
	{
		m_value = value;
	}

	void testInit()
	{
		Executor<TestSync> executor(this);

		async::as::Sync<Executor<TestSync>, Parameter, Parameter> async;
		async.setExecutor(executor);

		async.wait();
	}

	void testCall()
	{
		Executor<TestSync> executor(this);

		async::as::Sync<Executor<TestSync>, Parameter, Parameter> async;
		async.setExecutor(executor);

		TS_ASSERT_EQUALS(async.numBuffers(), 0);

		m_value = 0;

		async.wait();
		Parameter parameter;
		parameter.value = 42;
		async.call(parameter);

		async.wait(); // Make sure the call is finished
		TS_ASSERT_EQUALS(m_value, 42);

		parameter.value = 415;
		async.call(parameter);

		async.wait();
		TS_ASSERT_EQUALS(m_value, 415);
	}

	void testBuffer()
	{
		Executor<TestSync> executor(this);

		async::as::Sync<Executor<TestSync>, Parameter, Parameter> async;
		async.setExecutor(executor);

		int buffer = 2;
		TS_ASSERT_EQUALS(async.addBuffer(&buffer, sizeof(int)), 0);
		TS_ASSERT_EQUALS(async.numBuffers(), 1);

		TS_ASSERT_EQUALS(&buffer, async.buffer(0));

		async.wait();
	}

	void testRemoveBuffer()
	{
		Executor<TestSync> executor(this);

		async::as::Sync<Executor<TestSync>, Parameter, Parameter> async;
		async.setExecutor(executor);

		int buffer = 2;
		async.addBuffer(&buffer, sizeof(int));

		async.wait();

		async.removeBuffer(0);
		TS_ASSERT_EQUALS(static_cast<const void*>(0L), async.buffer(0));
	}

	void testManagedBuffer()
	{
		Executor<TestSync> executor(this);

		async::as::Sync<Executor<TestSync>, Parameter, Parameter> async;
		async.setExecutor(executor);

		async.addBuffer(0L, sizeof(int));

		*static_cast<int*>(async.managedBuffer(0)) = 2;

		TS_ASSERT_EQUALS(*static_cast<const int*>(async.buffer(0)), 2);

		async.wait();
	}
};
