/*! \file   DebugUtils.hpp
    \brief  compatibility of compilers
    \author Xavier Juvigny, ONERA
    \date   Jan. 12th 2005
*/

#ifndef _COMPILER_OPTIONCOMPILER_H
# define _COMPILER_OPTIONCOMPILER_H

// ========= Append a underscore or not for Fortran subroutines ========
#ifdef WIN32
# ifdef NB_NO_UNDERSCORE
#  define FORTRAN_DECL_WL(x_windows,x_linux) x_windows
# else
#   ifdef NB_DBLEUNDERSCORE
#     define FORTRAN_DECL_WL(x_windows,x_linux) x_windows##__
#   else
#     define FORTRAN_DECL_WL(x_windows,x_linux) x_windows##_
#   endif
# endif
#else
# ifdef NB_NO_UNDERSCORE
#  define FORTRAN_DECL_WL(x_windows,x_linux) x_linux
# else
#   ifdef NB_DBLEUNDERSCORE
#     define FORTRAN_DECL_WL(x_windows,x_linux) x_linux##__
#   else
#     define FORTRAN_DECL_WL(x_windows,x_linux) x_linux##_
#   endif
# endif
#endif

# ifdef NB_NO_UNDERSCORE
#  define FORTRAN_DECL(x) x
# else
#   ifdef NB_DBLEUNDERSCORE
#     define FORTRAN_DECL(x) x##__
#   else
#     define FORTRAN_DECL(x) x##_
#   endif
# endif

#ifdef _MSC_VER
#  ifdef _DLL
#    ifdef DISSECTION_EXPORTS // when building DLL
#      define DISSECTION_API __declspec(dllexport)
#    else // when client uses DLL
#      define DISSECTION_API __declspec(dllimport)
#    endif
#  else 
#    define DISSECTION_API
#  endif 
#else
#  define DISSECTION_API
#endif

#endif
