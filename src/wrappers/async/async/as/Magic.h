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

#ifndef ASYNC_AS_MAGIC_H
#define ASYNC_AS_MAGIC_H

namespace async
{

namespace as
{

/**
 * @class      : ASYNC_HAS_MEM_FUNC
 * @brief      : This macro will be used to check if a class has a particular
 * member function implemented in the public section or not.
 * @param func : Name of Member Function
 * @param name : Name of struct which is going to be run the test for
 * the given particular member function name specified in func
 * @param return_type: Return type of the member function
 * @param ellipsis(...) : Since this is macro should provide test case for every
 * possible member function we use variadic macros to cover all possibilities
 *
 * source: http://stackoverflow.com/questions/257288/is-it-possible-to-write-a-c-template-to-check-for-a-functions-existence
 * Modified to allow for an additional template function.
 */
#define ASYNC_HAS_MEM_FUNC_T1(func, name, T1, return_type, ...)   \
   template <typename T, typename T1>                             \
   struct name {                                                  \
      typedef return_type (T::*Sign)(__VA_ARGS__);                \
      typedef char yes[1];                                        \
      typedef char no[2];                                         \
      template <typename U, U>                                    \
      struct type_check;                                          \
      template <typename _1>                                      \
      static yes& chk(type_check<Sign, &_1::func>*);              \
      template <typename>                                         \
      static no& chk(...);                                        \
      static bool const value = sizeof(chk<T>(0)) == sizeof(yes); \
   }

template<bool C, typename T = void>
struct enable_if {
	typedef T type;
};

template<typename T>
struct enable_if<false, T> { };

}

}

#endif // ASYNC_AS_MAGIC_H