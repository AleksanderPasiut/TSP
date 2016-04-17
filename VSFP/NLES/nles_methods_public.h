#pragma once

#include <stdexcept>

template <typename T>
T NLES<T>::bisection_solver(T(*function)(const T&), T lower_limit, T higher_limit, const T& precision, unsigned iterations_limit, unsigned* iterations_counter)
{
	unsigned counter = 0;
	T mid;
	if (function(lower_limit) < 0)
		while (higher_limit - lower_limit > precision)
		{
			mid = (higher_limit+lower_limit)/2;
			T mid_value = function(mid);

			if (mid_value > 0)
				higher_limit = mid;
			else lower_limit = mid;

			if (++counter == iterations_limit)
				throw std::logic_error("Iterations limit reached.");
		}
	else
		while (higher_limit - lower_limit > precision)
		{
			mid = (higher_limit+lower_limit)/2;
			T mid_value = function(mid);

			if (mid_value < 0)
				higher_limit = mid;
			else lower_limit = mid;

			if (++counter == iterations_limit)
				throw std::logic_error("Iterations limit reached.");
		}

	if (iterations_counter)
		*iterations_counter = counter;

	return (higher_limit+lower_limit)/2;
}

template <typename T>
T NLES<T>::tangent_solver(T(*function)(const T&), T(*derivative)(const T&), T approx, const T& precision, unsigned iterations_limit, unsigned* iterations_counter)
{
	unsigned counter = 0;
	for(;;)
	{
		T tmp = function(approx)/derivative(approx);

		if (abs(tmp) < precision)
			break;

		approx -= tmp;

		if (++counter == iterations_limit)
			throw std::logic_error("Iterations limit reached.");
	}

	if (iterations_counter)
		*iterations_counter = counter;

	return approx;
}

template <typename T>
T NLES<T>::regula_falsi_solver(T(*function)(const T&), T lower_limit, T higher_limit, const T& precision, unsigned iterations_limit, unsigned* iterations_counter)
{
	unsigned counter = 0;
	const T f_lower_limit = function(lower_limit);
	const T f_higher_limit = function(higher_limit);
	T approx = (lower_limit*f_higher_limit-higher_limit*f_lower_limit)/(f_higher_limit-f_lower_limit);

	for(;;)
	{
		T f_approx = function(approx);

		if (abs(f_approx) < precision)
			break;

		if (same_sign(f_lower_limit, f_approx))
			approx -= f_approx/(f_higher_limit-f_approx) * (higher_limit-approx);
		else approx -= f_approx/(f_lower_limit-f_approx) * (lower_limit-approx);

		if (++counter == iterations_limit)
			throw std::logic_error("Iterations limit reached.");
	}

	if (iterations_counter)
		*iterations_counter = counter;

	return approx;
}

template <typename T>
T NLES<T>::secant_solver(T(*function)(const T&), T lower_limit, T higher_limit, const T& precision, unsigned iterations_limit, unsigned* iterations_counter)
{
	unsigned counter = 0;
	T* tvb = &lower_limit; // two values back
	T* ovb = &higher_limit; // one value back
	for(;;)
	{
		T f_ovb = function(*ovb);

		if (abs(f_ovb) < precision)
			break;

		T f_tvb = function(*tvb);
		*tvb = *ovb - (f_ovb*(*ovb-*tvb))/(f_ovb-f_tvb);

		T* tmp = tvb;
		tvb = ovb;
		ovb = tmp;

		if (++counter == iterations_limit)
			throw std::logic_error("Iterations limit reached.");
	}

	if (iterations_counter)
		*iterations_counter = counter;

	return *ovb;
}
