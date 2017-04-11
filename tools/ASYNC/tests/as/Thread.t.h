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

#include <pthread.h>
#include <cstdlib>
#include <ctime>
#include <sys/sysinfo.h>

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
		// Lock the variable to test multiple threads at once
		pthread_spin_lock(&m_lock);
		m_value += value;
		pthread_spin_unlock(&m_lock);
	}

	void setUp()
	{
		pthread_spin_init(&m_lock, PTHREAD_PROCESS_PRIVATE);
		m_value = 0;
	}

	void testInit()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		Parameter param;
		param.value = 32;
		async.callInit(param);

		TS_ASSERT_EQUALS(m_value, 32);

		async.wait();
	}

	void testInitBuffer()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		int buffer;
		async.addSyncBuffer(&buffer, sizeof(int));

		TS_ASSERT_EQUALS(async.numBuffers(), 1);
		TS_ASSERT_EQUALS(&buffer, async.buffer(0));
		TS_ASSERT_EQUALS(async.bufferSize(0), sizeof(int));

		async.wait();
	}

	void testCall()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		TS_ASSERT_EQUALS(async.numBuffers(), 0);

		async.wait();
		Parameter parameter;
		parameter.value = 42;
		async.call(parameter);

		async.wait(); // Make sure the call is finished
		TS_ASSERT_EQUALS(m_value, 42);

		parameter.value = 415;
		async.call(parameter);

		async.wait();
		TS_ASSERT_EQUALS(m_value, 42+415);
	}

	void testAffinity()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		// Set affinity (choose a random CPU)
		srand(time(0L));
		int cpu = rand() % get_nprocs();
		cpu_set_t cpuSet;
		CPU_ZERO(&cpuSet);
		CPU_SET(cpu, &cpuSet);
		async.setAffinity(cpuSet);

		async.wait();
		Parameter parameter;
		async.call(parameter);

		async.wait();

		TS_ASSERT_EQUALS(executor.cpu(), cpu);
	}

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
	}

	void testBuffer2()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		int buffer0 = 42;
		TS_ASSERT_EQUALS(async.addBuffer(&buffer0, sizeof(int)), 0);
		int buffer1 = 12;
		TS_ASSERT_EQUALS(async.addBuffer(&buffer1, sizeof(int)), 1);

		async.wait();
		async.sendBuffer(0, sizeof(int));
		TS_ASSERT_EQUALS(*reinterpret_cast<const int*>(async.buffer(0)), 42);

		async.sendBuffer(1, sizeof(int));
		TS_ASSERT_EQUALS(*reinterpret_cast<const int*>(async.buffer(1)), 12);
	}

	void testSyncBuffer()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		int buffer = 42;
		async.addSyncBuffer(&buffer, sizeof(int));

		async.wait();

		async.sendBuffer(0, sizeof(int));
		TS_ASSERT_EQUALS(async.buffer(0), &buffer);
	}

	void testManagedBuffer()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		async.addBuffer(0L, sizeof(int));

		TS_ASSERT_DIFFERS(async.managedBuffer(0), static_cast<void*>(0L));

		async.wait();

		*static_cast<int*>(async.managedBuffer(0)) = 32;
		async.sendBuffer(0, sizeof(int));

		TS_ASSERT_EQUALS(*reinterpret_cast<const int*>(async.buffer(0)), 32);
	}

	void testRemoveBuffer()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async;
		async.setExecutor(executor);

		int buffer0 = 42;
		async.addSyncBuffer(&buffer0, sizeof(int));
		int buffer1 = 12;
		async.addBuffer(&buffer1, sizeof(int));

		async.sendBuffer(0, sizeof(int));

		Parameter parameter;
		async.callInit(parameter);

		async.wait();

		async.removeBuffer(0);
		TS_ASSERT_EQUALS(async.buffer(0), static_cast<const void*>(0L));
		logInfo() << async.buffer(1);
		TS_ASSERT_DIFFERS(async.buffer(1), static_cast<const void*>(0L));

		async.sendBuffer(1, sizeof(int));

		async.call(parameter);

		async.wait();

		logInfo() << async.buffer(1);

		async.removeBuffer(1);
		TS_ASSERT_EQUALS(async.buffer(1), static_cast<const void*>(0L));
	}

	void testMultiple()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async1;
		async1.setExecutor(executor);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async2;
		async2.setExecutor(executor);

		async1.wait();
		async2.wait();

		Parameter parameter;
		parameter.value = 42;
		async1.call(parameter);

		parameter.value = 13;
		async2.call(parameter);

		async1.wait();
		async2.wait();

		TS_ASSERT_EQUALS(m_value, 42+13);
	}

	void testMultipleBuffer()
	{
		Executor<TestThread> executor(this);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async1;
		async1.setExecutor(executor);

		int buffer = 41;
		TS_ASSERT_EQUALS(async1.addBuffer(&buffer, sizeof(int)), 0);

		async::as::Thread<Executor<TestThread>, Parameter, Parameter> async2;
		async2.setExecutor(executor);

		TS_ASSERT_EQUALS(async2.addBuffer(&buffer, sizeof(int)), 0);

		async1.wait();
		async2.wait();

		async1.sendBuffer(0, sizeof(int));
		TS_ASSERT_EQUALS(*reinterpret_cast<const int*>(async1.buffer(0)), 41);

		Parameter parameter;
		parameter.value = 42;
		async1.call(parameter);

		parameter.value = 13;
		async2.call(parameter);

		async1.wait();
		async2.wait();

		TS_ASSERT_EQUALS(m_value, 42+13);
	}
};
