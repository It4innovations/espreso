/**
 * @file
 *  This file is part of UTILS
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2014-2015, Technische Universitaet Muenchen.
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

#ifndef UTILS_PATH_H
#define UTILS_PATH_H

#include <string>
#include <sys/stat.h>

#include "utils/stringutils.h"

namespace utils
{

/**
 * Manipulates file/directory names and paths
 */
class Path
{
private:
	std::string m_path;

public:
	Path()
	{
	}

	Path(const char* path)
		: m_path(path)
	{
		init();
	}

	Path(const std::string &path)
		: m_path(path)
	{
		init();
	}

	/**
	 * @return The string representing the current path
	 */
	operator std::string() const
	{
		return m_path;
	}

	/**
	 * @return The basename of the path
	 */
	std::string basename() const
	{
		const size_t lastSlash = m_path.find_last_of(separators());
		if (lastSlash == std::string::npos)
			return m_path;

		return m_path.substr(lastSlash+1);
	}

	/**
	 * @return The directory name of the path
	 */
	std::string dirname() const
	{
		const size_t lastSlash =  m_path.find_last_of(separators());
		if (lastSlash == std::string::npos)
			return "";

		return m_path.substr(0, lastSlash);
	}

	/**
	 * @return The directory of the path
	 */
	Path dir() const
	{
		return Path(dirname());
	}

	/**
	 * @return True of the path is a file/directory
	 */
	bool exists() const
	{
		struct stat buffer;
		return (stat(m_path.c_str(), &buffer) == 0);
	}

	/**
	 * Joins two paths
	 *
	 * @param other The path that should be appended (has to be relative)
	 * @return A path where other is appended to the current path
	 */
	Path operator+(const Path& other) const
	{
		if (m_path.empty())
			return other;
		if (other.m_path.empty())
			return *this;

		return Path(m_path + SEPARATOR + other.m_path);
	}

public:
	static const char* separator()
	{
		static const std::string sep(1, SEPARATOR);
		return sep.c_str();
	}

	/**
	 * @return A string containing all possible separators
	 */
	static const char* separators()
	{
#if __unix__
		return "/";
#else // __unix__
		return "\\/";
#endif // __unix__
	}

private:
	void init()
	{
		// Remove trailing separator
		if (StringUtils::endsWith(m_path, separator()))
			StringUtils::replaceLast(m_path, separator(), "");
	}

public:
#ifdef __unix__
	static const char SEPARATOR = '/';
#else // __unix__
	static const char SEPARATOR = '\\';
#endif // __unix__

};

}

#endif // UTILS_PATH_H
