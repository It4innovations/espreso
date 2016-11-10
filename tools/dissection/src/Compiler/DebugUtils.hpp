/*! \file   DebugUtils.hpp
    \brief  Some macros and functions to help debug session
    \author Xavier Juvigny, ONERA
    \date   Jan. 19th 2005
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _DISSECTION_COMPILER_DEBUGUTILS_HPP_
#define _DISSECTION_COMPILER_DEBUGUTILS_HPP_
#include <cassert>

#if defined(DISSECTION_DEBUG)
#  define CHECK(o,msg)				\
    assert((o)&&(msg))
#else
#  define CHECK(o,msg) 
#endif

#if defined(DISSECTION_DEBUG)				
#  define DBG_PRINT printf
#else
#  define DBG_PRINT //
#endif

#if defined(TRACE)
#  define TRACE printf
#else
#  define TRACE //
#endif

#endif
