/**
 * @file
 *  This file is part of UTILS
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2013-2016, Technische Universitaet Muenchen.
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

#ifndef UTILS_STRINGUTILS_H
#define UTILS_STRINGUTILS_H

#include <algorithm>
#include <functional>
#include <cctype>
#include <cstring>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

/**
 * A collection of useful utility functions
 */
namespace utils
{

/**
 * A collection of useful string functions based on std::string
 */
class StringUtils
{
public:

	/**
	 * Replaces first occurrence of from in str with to
	 *
	 * Taken from http://stackoverflow.com/questions/3418231/c-replace-part-of-a-string-with-another-string
	 */
	static bool replace(std::string& str, const std::string& from, const std::string& to) {
		size_t start_pos = str.find(from);
		if(start_pos == std::string::npos)
			return false;
		str.replace(start_pos, from.length(), to);
		return true;
	}

	/**
	 * Replaces last occurrence of from in str with to
	 */
	static bool replaceLast(std::string& str, const std::string& from, const std::string& to) {
		size_t start_pos = str.rfind(from);
		if (start_pos == std::string::npos)
			return false;
		str.replace(start_pos, from.length(), to);
		return true;
	}

	static bool startsWith(const std::string &str, const std::string &prefix)
	{
		return str.size() >= prefix.size() &&
				str.compare(0, prefix.size(), prefix) == 0;
	}

	/**
	 * Checks a string for specific suffix
	 *
	 * Taken from: http://stackoverflow.com/questions/20446201/how-to-check-if-string-ends-with-txt
	 */
	static bool endsWith(const std::string &str, const std::string &suffix)
	{
		return str.size() >= suffix.size() &&
				str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
	}

	/**
	 * Converts arbitrary datatypes (all datatypes which support the << stream
	 * operator) into std::string
	 */
	template<typename T>
	static std::string toString(T value)
	{
		std::ostringstream ss;
		ss << value;
		return ss.str();
	}

	/**
	 * Converts strings to arbitrary datatypes (using the << stream operator)
	 *
	 * @param str The string that should be converted
	 */
	template<typename T>
	static T parse(const std::string &str)
	{
		T result;
		std::istringstream(str) >> result;

		return result;
	}

	/**
	 * Converts strings to arbitrary datatypes
	 *
	 * @param str The string that should be converted
	 * @param advanced True of advanced conversions should be enabled
	 */
	template<typename T>
	static T parse(const std::string &str, bool adavanced)
	{
		// By default the advanced mode is disabled for all datatypes
		return parse<T>(str);
	}

	template<typename T> inline
	static std::vector<T> parseArray(const std::string &str)
	{
		std::vector<T> elems;
		std::istringstream f(str);
		std::string s;
		while (std::getline(f, s, ':'))
			elems.push_back(parse<T>(s));
		
		return elems;
	}

	/**
	 * Converts null terminated string to upper case
	 */
	static void toUpper(char* s)
	{
		for (size_t i = 0; s[i]; i++)
			s[i] = toupper(s[i]);
	}

	/**
	 * Converts std::string to upper case
	 */
	static void toUpper(std::string &str)
	{
		for (size_t i = 0;  str[i]; i++)
			str[i] = toupper(str[i]);
	}

	/**
	 * Converts null terminated string to lower case
	 */
	static void toLower(char* s)
	{
		for (size_t i = 0; s[i]; i++)
			s[i] = tolower(s[i]);
	}

	/**
	 * Converts std::string to lower case
	 */
	static void toLower(std::string &str)
	{
		for (size_t i = 0; str[i]; i++)
			str[i] = tolower(str[i]);
	}

	/**
	 * Trims from start
	 *
	 * Taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
	 */
	static std::string &ltrim(std::string &s) {
	        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	        return s;
	}

	/**
	 * Trims from end
	 *
	 * Taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
	 */
	static std::string &rtrim(std::string &s) {
	        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	        return s;
	}

	/**
	 * Trims from both ends
	 *
	 * Taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
	 */
	static std::string &trim(std::string &s) {
	        return ltrim(rtrim(s));
	}

	/**
	 * Join vector
	 *
	 * Taken from http://dracoater.blogspot.com/2010/03/joining-vector-and-splitting-string-in.html
	 */
	template<typename T>
	static std::string join(const std::vector<T> &v, const std::string &token) {
		std::ostringstream result;
		for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); i++){
			if (i != v.begin())
				result << token;
			result << *i;
		}

		return result.str();
	}

	/**
	 * Split a string
	 *
	 * Taken from http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
	 */
	static std::vector<std::string> split(const std::string &str, char delim)
	{
		std::vector<std::string> elems;
		std::stringstream ss(str);

		std::string item;
		while (getline(ss, item, delim))
			elems.push_back(item);

		return elems;
	}
};

template<> inline
std::string StringUtils::parse(const std::string &str)
{
	return str;
}

template<> inline
const char* StringUtils::parse(const std::string &str)
{
	return str.c_str();
}

template<> inline
bool StringUtils::parse(const std::string &str, bool advanced)
{
	std::string s = str;
	toLower(s);

	if (s == "on"
			|| s == "yes")
		return true;

	return parse<bool>(str);
}

}

#endif /* UTILS_STRINGUTILS_H */
