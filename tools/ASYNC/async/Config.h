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

#ifndef ASYNC_CONFIG_H
#define ASYNC_CONFIG_H

#include <string>

#include "utils/env.h"
#include "utils/logger.h"
#include "utils/stringutils.h"

namespace async
{

/**
 * The asynchnchronous mode that should be used
 */
enum Mode
{
	SYNC,
	THREAD,
	MPI
};

/**
 * @warning Overwriting the default values from the environment is only allowed
 *  before any other ASYNC class is created. Otherwise the behavior is undefined.
 */
class Config
{
private:
	Mode m_mode;
	
	int m_pinCore;
	
	unsigned int m_groupSize;
	
	bool m_asyncCopy;
	
	size_t m_alignment;
	
private:
	Config()
		: m_mode(str2mode(utils::Env::get<const char*>("ASYNC_MODE", "SYNC"))),
		m_pinCore(utils::Env::get<int>("ASYNC_PIN_CORE", -1)),
		m_groupSize(m_mode == MPI ? utils::Env::get("ASYNC_GROUP_SIZE", 64) : 1),
		m_asyncCopy(utils::Env::get<int>("ASYNC_MPI_COPY", 0)),
		m_alignment(utils::Env::get<size_t>("ASYNC_BUFFER_ALIGNMENT", 0))
	{ }
	
public:
	static Mode mode()
	{
		return instance().m_mode;
	}

	static int getPinCore()
	{
		return instance().m_pinCore;
	}

	static unsigned int groupSize()
	{
		return instance().m_groupSize;
	}

	static bool useAsyncCopy()
	{
		return instance().m_asyncCopy;
	}

	static size_t alignment()
	{
		return instance().m_alignment;
	}

	static size_t maxSend()
	{
		return utils::Env::get<size_t>("ASYNC_MPI_MAX_SEND", 1UL<<30);
	}
	
	static void setMode(Mode mode)
	{
		instance().m_mode = mode;
	}
	
	static void setPinCore(int pinCore)
	{
		instance().m_pinCore = pinCore;
	}
	
	static void setUseAsyncCopy(bool useAsyncCopy)
	{
		instance().m_asyncCopy = useAsyncCopy;
	}
	
	static void setGroupSize(unsigned int groupSize)
	{
		if (instance().mode() == MPI)
			instance().m_groupSize = groupSize;
	}
	
	static void setAlignment(size_t alignment)
	{
		instance().m_alignment = alignment;
	}

private:
	static Config& instance()
	{
		static Config config;
		return config;
	}
	
	static Mode str2mode(const char* mode)
	{
		std::string strMode(mode);
		utils::StringUtils::toUpper(strMode);

		if (strMode == "THREAD")
			return THREAD;
		if (strMode == "MPI") {
#ifdef USE_MPI
			return MPI;
#else // USE_MPI
			logError() << "Asynchronous MPI is not supported without MPI";
#endif // USE_MPI
		}
		if (strMode != "SYNC")
			logWarning() << "Unknown mode" << utils::nospace << strMode << ". Using synchronous mode.";
		return SYNC;
	}
};

}

#endif // ASYNC_CONFIG_H