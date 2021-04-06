/**
 * @file
 *  This file is part of UTILS
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2014-2016, Technische Universitaet Muenchen.
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

#include "utils/stringutils.h"

using namespace utils;

class TestStringUtils : public CxxTest::TestSuite
{
public:
	void testStartsWith()
	{
		TS_ASSERT(StringUtils::startsWith("abcde", "abc"));
		TS_ASSERT(!StringUtils::startsWith("abcde", "abd"));
	}

	void testEndsWith()
	{
		TS_ASSERT(StringUtils::endsWith("abcde", "cde"));
		TS_ASSERT(!StringUtils::endsWith("abcde", "abc"));
	}

	void testParse()
	{
		// Normal parser
		// TODO more tests
		TS_ASSERT_EQUALS(StringUtils::parse<int>("-1", 0), -1);

		// Advanced parser
		TS_ASSERT(StringUtils::parse<bool>("on", true));
		TS_ASSERT(StringUtils::parse<bool>("yes", true));
		TS_ASSERT(StringUtils::parse<bool>("on", true));
		TS_ASSERT(!StringUtils::parse<bool>("off", true));
		TS_ASSERT(!StringUtils::parse<bool>("abc", true));
	}
	
	void testParseArray()
	{
		// TODO more tests
		std::vector<int> result = StringUtils::parseArray<int>("1:2:3");
		
		TS_ASSERT_EQUALS(result[0], 1);
		TS_ASSERT_EQUALS(result[1], 2);
		TS_ASSERT_EQUALS(result[2], 3);
	}
};
