// This file provides MATRIX_SoLE template structure which extends MATRIX template class
// with functions to perform more advanced mathematical operations such as determinant
// computation, SoLE (System of Linear Equations) solvers, etc.

#pragma once

#include "..\matrix.h"

template<typename T>
class MATRIX_SoLE
{
protected:
	static MATRIX<T>& diagonal(MATRIX<T>&);
	static MATRIX<T>& nondiagonal(MATRIX<T>&);
	static MATRIX<T>& upper(MATRIX<T>&);
	static MATRIX<T>& nonupper(MATRIX<T>&);
	static MATRIX<T>& row_echelon_form(MATRIX<T>&);
	static MATRIX<T>& jordan_elimination(MATRIX<T>&);
	static MATRIX<T>& diagonal_matrix_inverse(MATRIX<T>&);
	static MATRIX<T> lower_triangular_matrix_inverse(const MATRIX<T>&);
	static MATRIX<T>& shuffle_for_diagonal_domination(MATRIX<T>&, MATRIX<T>&);

	static T diagonal_elements_product(const MATRIX<T>&);

public:
	static MATRIX<T> matrix_inverse(const MATRIX<T>& matrix);
	static MATRIX<T> gauss_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector);
	static MATRIX<T> jordan_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector);
	static MATRIX<T> jacoby_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector, T precision);
	static MATRIX<T> gauss_seidl_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector, T precision);

	static T compute_determinant(const MATRIX<T>& matrix);
	static T max_norm(const MATRIX<T>&);
};

#include "matrix_SoLE_methods_protected.h"
#include "matrix_SoLE_methods_public.h"
