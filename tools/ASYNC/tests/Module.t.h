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

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <cxxtest/TestSuite.h>

#include <sched.h>

#include "async/Config.h"
#include "async/Dispatcher.h"
#include "async/Module.h"

struct Param
{
	int step;
};

class Module : private async::Module<Module, Param, Param>
{
public:
	bool m_setUp;
	bool m_execInit;
	bool m_exec;
	bool m_tearDown;

	int m_cpu;

public:
	Module()
		: m_setUp(false),
		  m_execInit(false),
		  m_exec(false),
		  m_tearDown(false)
	{
	}

	void run()
	{
		init();

		Param param;
		callInit(param);

		wait();

		call(param);

		wait();

		finalize();
	}

	void norun()
	{
		init();
		finalize();
	}

	void execInit(const async::ExecInfo &info, const Param &param)
	{
		TS_ASSERT_EQUALS(info.numBuffers(), 0);
		m_execInit = true;
	}

	void exec(const Param &param)
	{
		m_exec = true;

		m_cpu = sched_getcpu();
	}

	void setUp()
	{
		setExecutor(*this);
		m_setUp = true;
	}

	void tearDown()
	{
		m_tearDown = true;
	}
};

class BufferModule : private async::Module<BufferModule, Param, Param>
{
public:
	unsigned int m_initBufferSize;
	unsigned int m_bufferSize;
	unsigned int m_cloneBufferSize;
	unsigned int m_cloneSyncBufferSize;
	int m_buffer;
	int m_managedBuffer;

public:
	BufferModule()
		: m_initBufferSize(0),
		  m_bufferSize(0),
		  m_cloneBufferSize(0),
		  m_cloneSyncBufferSize(0)
	{
	}

	void run()
	{
		init();

		int initBuffer[2];
		addBuffer(initBuffer, 2*sizeof(int));

		int buffer = 42;
		addBuffer(&buffer, sizeof(int));

		int cloneBuffer = 1;
		addBuffer(&cloneBuffer, sizeof(int), true);

		long addCloneSyncBuffer = 2;
		addSyncBuffer(&addCloneSyncBuffer, sizeof(long), true);

		addBuffer(0L, 2*sizeof(int));

		int buffer2 = 43;
		addBuffer(&buffer2, sizeof(int));

		int* managedBuffer = async::Module<BufferModule, Param, Param>::managedBuffer<int*>(4);
		TS_ASSERT_DIFFERS(managedBuffer, static_cast<int*>(0L));

		Param param;
		callInit(param);

		wait();

		sendBuffer(1, sizeof(int));

		managedBuffer[0] = 5;
		managedBuffer[1] = 5;

		sendBuffer(4, 2*sizeof(int));

		sendBuffer(5);

		call(param);

		if (async::Config::mode() == async::MPI) {
			// Set the params on non-executors
			execInit(param);
			_exec(param, false);
		}

		wait();

		finalize();
	}

	void execInit(const Param &param)
	{
		m_initBufferSize = bufferSize(0);
	}

	void exec(const async::ExecInfo &info, const Param &param)
	{
		_exec(param);
	}

	void setUp()
	{
		setExecutor(*this);
	}

	void tearDown()
	{
	}

private:
	void _exec(const Param &param, bool isExecutor = true)
	{
		m_bufferSize = bufferSize(1);
		for (unsigned int i = 0; i < m_bufferSize/sizeof(int); i++) {
			TS_ASSERT_EQUALS(42, *(static_cast<const int*>(buffer(1))+i));
		}
		m_cloneBufferSize = bufferSize(2);
		m_cloneSyncBufferSize = bufferSize(3);

		unsigned int bufferSize2 = bufferSize(5);
		TS_ASSERT_EQUALS(bufferSize2, m_bufferSize);
		for (unsigned int i = 0; i < bufferSize2/sizeof(int); i++) {
			TS_ASSERT_EQUALS(43, *(static_cast<const int*>(buffer(5))+i));
		}

		if (isExecutor) {
			unsigned int managedBufferSize = bufferSize(4);
			for (unsigned int i = 0; i < managedBufferSize/sizeof(int); i++) {
				TS_ASSERT_EQUALS(*(static_cast<const int*>(buffer(4))+i), 5);
			}
		}
	}
};



class ResizeBufferModule : private async::Module<ResizeBufferModule, Param, Param>
{
public:
	unsigned int m_buffer0Size[2];
	unsigned int m_buffer1Size[2];

public:
	ResizeBufferModule()
	{
	}

	void run()
	{
		init();

		int buffer0 = 42;
		addBuffer(&buffer0, sizeof(int));

		int buffer1[2] = {1, 2};
		addBuffer(buffer1, 2*sizeof(int), true);

		Param param;
		callInit(param);

		wait();

		sendBuffer(0, sizeof(int));
		sendBuffer(1, 2*sizeof(int));

		param.step = 0;
		call(param);

		if (async::Config::mode() == async::MPI) {
			// Set the params on non-executors
			_exec(param, false);
		}

		wait();
		
		int buffer2[2] = {4, 3};
		resizeBuffer(0, buffer2, 2*sizeof(int));
		
		resizeBuffer(1, buffer1, sizeof(int));
		
		sendBuffer(0);
		sendBuffer(1);
		
		param.step = 1;
		call(param);

		if (async::Config::mode() == async::MPI) {
			// Set the params on non-executors
			_exec(param, false);
		}
		
		wait();

		finalize();
	}

	void execInit(const Param &param)
	{
	}

	void exec(const async::ExecInfo &info, const Param &param)
	{
		_exec(param);
	}

	void setUp()
	{
		setExecutor(*this);
	}

	void tearDown()
	{
	}

private:
	void _exec(const Param &param, bool isExecutor = true)
	{
		m_buffer0Size[param.step] = bufferSize(0);
		m_buffer1Size[param.step] = bufferSize(1);
		switch (param.step) {
		case 0:
			TS_ASSERT_EQUALS(*static_cast<const int*>(buffer(0)), 42);
			TS_ASSERT_EQUALS(*static_cast<const int*>(buffer(1)), 1);
			TS_ASSERT_EQUALS(*(static_cast<const int*>(buffer(1))+1), 2);
			break;
		case 1:
			TS_ASSERT_EQUALS(*static_cast<const int*>(buffer(0)), 4);
			TS_ASSERT_EQUALS(*(static_cast<const int*>(buffer(0))+1), 3);
			TS_ASSERT_EQUALS(*static_cast<const int*>(buffer(1)), 1);
			break;
		}
	}
};

class RemoveBufferModule : private async::Module<RemoveBufferModule, Param, Param>
{
public:
	unsigned int m_buffer0Size;
	unsigned int m_buffer1Size;
	unsigned int m_managedBufferSize;

public:
	RemoveBufferModule()
		: m_buffer0Size(42),
		  m_buffer1Size(42)
	{
	}

	void run()
	{
		init();

		int buffer0 = 3;
		addBuffer(&buffer0, sizeof(int));

		int buffer1 = 42;
		addBuffer(&buffer1, sizeof(int));

		addBuffer(0L, sizeof(int));

		Param param;
		callInit(param);

		wait();

		removeBuffer(0);
		removeBuffer(2);
		sendBuffer(1, sizeof(int));

		call(param);

		if (async::Config::mode() == async::MPI) {
			// Set the params on non-executors
			execInit(param);
			exec(param);
		}

		wait();

		finalize();
	}

	void execInit(const Param &param)
	{
	}

	void exec(const Param &param)
	{
		m_buffer0Size = bufferSize(0);
		m_buffer1Size = bufferSize(1);
		m_managedBufferSize = bufferSize(2);
	}

	void setUp()
	{
		setExecutor(*this);
	}

	void tearDown()
	{
	}
};

/**
 * Test for {@link Dispatcher} and {@link Module} since
 * they work closely together.
 */
class TestModule : public CxxTest::TestSuite
{
	int m_rank;
	int m_size;

public:
	void setUp()
	{
#ifdef USE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
		MPI_Comm_size(MPI_COMM_WORLD, &m_size);
#endif // USE_MPI
	}

	void testModuleDispatcher()
	{
		async::Dispatcher dispatcher;
#ifdef USE_MPI
		dispatcher.setCommunicator(MPI_COMM_WORLD);
#endif // USE_MPI

		unsigned int groupSize = dispatcher.groupSize();
		if (async::Config::mode() == async::MPI) {
			TS_ASSERT_EQUALS(groupSize, 64); // the default
		} else {
			TS_ASSERT_EQUALS(groupSize, 1);
		}

		Module module;

		dispatcher.init();

		if (dispatcher.dispatch()) {
			TS_ASSERT(!dispatcher.isExecutor());

			module.run();
		} else {
			TS_ASSERT(dispatcher.isExecutor());

			if (async::Config::mode() == async::MPI) {
				TS_ASSERT_EQUALS(m_rank+1, m_size);
			} else {
				TS_FAIL("No executors in SYNC and THREAD mode!");
			}
		}

		dispatcher.finalize();

		TS_ASSERT(module.m_setUp);
		if (dispatcher.isExecutor()) {
			TS_ASSERT(module.m_execInit);
			TS_ASSERT(module.m_exec);
		}
		TS_ASSERT(module.m_tearDown);

		if (async::Config::mode() == async::THREAD) {
			TS_ASSERT_EQUALS(module.m_cpu, get_nprocs()-1);
		}
	}

	void testBuffer()
	{
		async::Dispatcher dispatcher;

		BufferModule module;

		dispatcher.init();

		if (dispatcher.dispatch())
			module.run();

		unsigned int initBufferSize = 2*sizeof(int);
		unsigned int bufferSize = sizeof(int);
		unsigned int cloneBufferSize = sizeof(int);
		unsigned int cloneSyncBufferSize = sizeof(long);

		if (dispatcher.isExecutor() && async::Config::mode() == async::MPI) {
			initBufferSize *= m_size-1;
			bufferSize *= m_size-1;
		}

		dispatcher.finalize();

		TS_ASSERT_EQUALS(module.m_initBufferSize, initBufferSize);
		TS_ASSERT_EQUALS(module.m_bufferSize, bufferSize);
		TS_ASSERT_EQUALS(module.m_cloneBufferSize, cloneBufferSize);
		TS_ASSERT_EQUALS(module.m_cloneSyncBufferSize, cloneSyncBufferSize);
	}

	void testResizeBuffer()
	{
		async::Dispatcher dispatcher;

		ResizeBufferModule module;

		dispatcher.init();

		if (dispatcher.dispatch())
			module.run();

		unsigned int bufferSize = sizeof(int);

		if (dispatcher.isExecutor() && async::Config::mode() == async::MPI) {
			bufferSize *= m_size-1;
		}

		dispatcher.finalize();

		TS_ASSERT_EQUALS(module.m_buffer0Size[0], bufferSize);
		TS_ASSERT_EQUALS(module.m_buffer0Size[1], 2*bufferSize);
		TS_ASSERT_EQUALS(module.m_buffer1Size[0], 2*sizeof(int));
		TS_ASSERT_EQUALS(module.m_buffer1Size[1], sizeof(int));
	}

	void testRemoveBuffer()
	{
		async::Dispatcher dispatcher;

		RemoveBufferModule module;

		dispatcher.init();

		if (dispatcher.dispatch())
			module.run();

		unsigned int buffer1Size = sizeof(int);

		if (dispatcher.isExecutor() && async::Config::mode() == async::MPI) {
			buffer1Size *= m_size-1;
		}

		dispatcher.finalize();

		TS_ASSERT_EQUALS(module.m_buffer0Size, 0);
		TS_ASSERT_EQUALS(module.m_buffer1Size, buffer1Size);
		TS_ASSERT_EQUALS(module.m_managedBufferSize, 0);
	}

	void testModulesActiveDisabled()
	{
		async::Dispatcher dispatcher;

		unsigned int groupSize = dispatcher.groupSize();
		if (async::Config::mode() == async::MPI) {
			TS_ASSERT_EQUALS(groupSize, 64); // the default
		} else {
			TS_ASSERT_EQUALS(groupSize, 1);
		}

		Module module0;
		Module module1;

		dispatcher.init();

		if (dispatcher.dispatch()) {
			// Run only one module
			module0.run();
		}

		dispatcher.finalize();
	}
};
