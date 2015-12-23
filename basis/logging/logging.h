
#ifndef BASIS_LOGGING_LOGGING_H_
#define BASIS_LOGGING_LOGGING_H_

namespace eslog {

struct Info1 {
	template<typename T>
	Info1& operator<<(const T& value)
	{
#if VERBOSE > 0
		std::cout << value;
#endif
		return *this;
	}

};

struct Info2 {
	template<typename T>
	Info2& operator<<(const T& value)
	{
#if VERBOSE > 1
		std::cout << value;
#endif
		return *this;
	}
};

struct Info3 {
	template<typename T>
	Info3& operator<<(const T& value)
	{
#if VERBOSE > 2
		std::cout << value;
#endif
		return *this;
	}
};

struct Debug {
	template<typename T>
	Debug& operator<<(const T& value)
	{
#ifdef DEBUG
		std::cout << value;
#endif
		return *this;
	}
};

struct Error {
	template<typename T>
	Error& operator<<(const T& value)
	{
		std::cerr << value;
		return *this;
	}
};

extern Info1 info1;
extern Info2 info2;
extern Info3 info3;
extern Debug debug;
extern Error error;
}


#endif /* BASIS_LOGGING_LOGGING_H_ */
