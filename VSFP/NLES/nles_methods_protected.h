#pragma once

#include <stdexcept>

template <typename T>
inline bool NLES<T>::same_sign(const T& arg1, const T& arg2)
{
	if (arg1 > 0 && arg2 > 0)
		return true;
	else if (arg1 < 0 && arg2 < 0)
		return true;
	else return false;
}