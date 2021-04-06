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

#include <cxxtest/TestSuite.h>

#include "utils/path.h"

using namespace utils;

class TestPath : public CxxTest::TestSuite
{
public:
	void testToString()
	{
		TS_ASSERT_EQUALS(std::string(Path("foo/")), "foo");
	}

	void testBasename()
	{
		TS_ASSERT_EQUALS(Path("foo/bar").basename(), "bar");
		TS_ASSERT_EQUALS(Path("foo").basename(), "foo")
	}

	void testDirname()
	{
		TS_ASSERT_EQUALS(Path("foo/bar").dirname(), "foo");
		TS_ASSERT_EQUALS(Path("foo/foo/bar").dirname(), "foo/foo");
		TS_ASSERT_EQUALS(Path("foo").dirname(), "");
	}

	void testDir()
	{
		TS_ASSERT_EQUALS(std::string(Path("foo/bar").dir()), "foo");
	}

	void testExists()
	{
		TS_ASSERT(Path("/dev/null").exists());
		TS_ASSERT(!Path("/dev/asdfasdf").exists());
	}

	void testOperatorPlus()
	{
		TS_ASSERT_EQUALS(std::string(Path("foo")+Path("bar")), "foo/bar");
		TS_ASSERT_EQUALS(std::string(Path("foo/")+Path("bar")), "foo/bar");
		TS_ASSERT_EQUALS(std::string(Path("foo")+Path("")), "foo");
	}
};
