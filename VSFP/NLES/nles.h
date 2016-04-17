// This folder provides NLES template class which consists of functions providing various
// Non-Linear Equation Solvers.

#pragma once

template <typename T>
class NLES
{
protected:
	static bool same_sign(const T&, const T&);

public:
	static T bisection_solver(	T(*function)(const T&),
								T lower_limit,
								T higher_limit,
								const T& precision,
								unsigned iterations_limit = 100,
								unsigned* iterations_counter = 0);

	static T tangent_solver(	T(*function)(const T&),
								T(*derivative)(const T&),
								T starting_point,
								const T& precision,
								unsigned iterations_limit = 100,
								unsigned* iterations_counter = 0);

	static T regula_falsi_solver(	T(*function)(const T&),
									T lower_limit,
									T higher_limit,
									const T& precision,
									unsigned iterations_limit = 100,
									unsigned* iterations_counter = 0);

	static T secant_solver(	T(*function)(const T&),
							T lower_limit,
							T higher_limit,
							const T& precision,
							unsigned iterations_limit = 100,
							unsigned* iterations_counter = 0);
};

#include "nles_methods_protected.h"
#include "nles_methods_public.h"