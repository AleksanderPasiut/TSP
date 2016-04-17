#pragma once

#include <stdexcept>

template <typename T>
MATRIX<T> MATRIX_SoLE<T>::matrix_inverse(const MATRIX<T>& matrix)
{
	if (matrix.rows() != matrix.cols())
		throw std::logic_error("Matrix is not square.");

	MATRIX<T> tmp = matrix;
	MATRIX<T> ret(matrix.rows(), matrix.cols());
	ret.identity();

	for (unsigned i = 0; i+1 < tmp.cols(); i++)
	{
		if (tmp.field(i,i) == 0)
			for (unsigned j = i+1; j < tmp.rows(); j++)
				if (tmp.field(j,i) != 0)
					for(unsigned k = 0; k < tmp.cols(); k++)
					{	T t = tmp.field(i,k);
						tmp.field(i,k) = tmp.field(j,k);
						tmp.field(j,k) = t;
						t = ret.field(i,k);
						ret.field(i,k) = ret.field(j,k);
						ret.field(j,k) = t; }

		if (tmp.field(i,i) != 0)
			for (unsigned j = i+1; j < tmp.rows(); j++)
			{
				T multiplier = tmp.field(j,i)/tmp.field(i,i);

				for (unsigned k = i; k < tmp.cols(); k++)
					tmp.field(j,k) -= tmp.field(i,k)*multiplier;

				for (unsigned k = 0; k < tmp.cols(); k++)
					ret.field(j,k) -= ret.field(i,k)*multiplier;
			}
	}

	if (diagonal_elements_product(tmp) == 0)
		throw std::logic_error("Matrix is not invertible.");
	
	for (unsigned i = tmp.cols()-1; i > 0; i--)
		for (unsigned j = 0; j < i; j++)
		{
			T multiplier = tmp.field(j,i)/tmp.field(i,i);
			for (unsigned k = 0; k < ret.cols(); k++)
				ret.field(j,k) -= multiplier*ret.field(i,k);
		}

	for (unsigned i = 0; i < tmp.rows(); i++)
		for (unsigned j = 0; j < tmp.cols(); j++)
			ret.field(i,j) /= tmp.field(i,i);

	return ret;
}

template <typename T>
MATRIX<T> MATRIX_SoLE<T>::gauss_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector)
{
	if (vector.cols() != 1)
		throw std::logic_error("Vector has more than one column.");

	if (matrix.cols() != matrix.rows())
		throw std::logic_error("Matrix is not invertible.");

	if (matrix.rows() != vector.rows())
		throw std::logic_error("Vector and matrix have different amount of rows.");

	MATRIX<T> tmp = matrix;
	tmp.col_combine(vector);

	row_echelon_form(tmp);

	if (diagonal_elements_product(tmp) == 0)
		throw std::logic_error("Matrix is not invertible.");

	MATRIX<T> ret(vector.rows(), vector.cols());

	ret.field(vector.rows()-1, 0) = tmp.field(vector.rows()-1, vector.rows())/tmp.field(vector.rows()-1, vector.rows()-1);

	if (vector.rows() >= 2)
		for (unsigned i = vector.rows()-2; ; i--)
		{
			ret.field(i,0) = tmp.field(i,vector.rows());

			for (unsigned j = i+1; j < vector.rows(); j++)
				ret.field(i,0) -= tmp.field(i,j)*ret.field(j,0);

			ret.field(i,0) /= tmp.field(i,i);

			if (!i)
				break;
		}

	return ret;
}

template <typename T>
MATRIX<T> MATRIX_SoLE<T>::jordan_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector)
{
	if (vector.cols() != 1)
		throw std::logic_error("Vector has more than one column.");

	if (matrix.cols() != matrix.rows())
		throw std::logic_error("Matrix is not invertible.");

	if (matrix.rows() != vector.rows())
		throw std::logic_error("Vector and matrix have different amount of rows.");

	MATRIX<T> tmp = matrix;
	tmp.col_combine(vector);

	jordan_elimination(tmp);

	if (diagonal_elements_product(tmp) == 0)
		throw std::logic_error("Matrix is not invertible.");

	MATRIX<T> ret(vector.rows(), vector.cols());

	for (unsigned i = 0; i < vector.rows(); i++)
		ret.field(i,0) = tmp.field(i,matrix.cols())/tmp.field(i,i);

	return ret;
}

template<typename T>
MATRIX<T> MATRIX_SoLE<T>::jacoby_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector, T precision)
{
	if (vector.cols() != 1)
		throw std::logic_error("Vector has more than one column.");

	if (matrix.cols() != matrix.rows())
		throw std::logic_error("Matrix is not invertible.");

	if (matrix.rows() != vector.rows())
		throw std::logic_error("Vector and matrix have different amount of rows.");

	if (precision < 0)
		throw std::logic_error("Precision cannot be negative.");

	MATRIX<T> m = matrix;
	MATRIX<T> v2 = vector;
	MATRIX<T> vec = diagonal_matrix_inverse(diagonal(shuffle_for_diagonal_domination(m, v2)))*v2;
	MATRIX<T> LU = matrix;
	m *= nondiagonal(LU)*(-1);

	MATRIX<T> m1 = vec;
	MATRIX<T> m2;

	MATRIX<T>* approx1 = &m1;
	MATRIX<T>* approx2 = &m2;

	for (;;)
	{
		*approx2 = m*(*approx1)+vec;

		if (max_norm(*approx2-*approx1) < precision)
			break;

		MATRIX<T>* tmp = approx1;
		approx1 = approx2;
		approx2 = tmp;
	}
	return *approx2;
}

template<typename T>
MATRIX<T> MATRIX_SoLE<T>::gauss_seidl_solver(const MATRIX<T>& matrix, const MATRIX<T>& vector, T precision)
{
	if (vector.cols() != 1)
		throw std::logic_error("Vector has more than one column.");

	if (matrix.cols() != matrix.rows())
		throw std::logic_error("Matrix is not invertible.");

	if (matrix.rows() != vector.rows())
		throw std::logic_error("Vector and matrix have different amount of rows.");

	if (precision < 0)
		throw std::logic_error("Precision cannot be negative.");

	MATRIX<T> mGS = matrix;
	MATRIX<T> v2 = vector;
	MATRIX<T> vec = lower_triangular_matrix_inverse(nonupper(shuffle_for_diagonal_domination(mGS, v2)))*v2;
	MATRIX<T> up = matrix;
	mGS *= upper(up) * (-1);

	MATRIX<T> m1 = vec;
	MATRIX<T> m2;

	MATRIX<T>* approx1 = &m1;
	MATRIX<T>* approx2 = &m2;

	for (;;)
	{
		*approx2 = mGS*(*approx1)+vec;

		if (max_norm(*approx2-*approx1) < precision)
			break;

		MATRIX<T>* tmp = approx1;
		approx1 = approx2;
		approx2 = tmp;
	}
	return *approx2;
}

template <typename T>
T MATRIX_SoLE<T>::compute_determinant(const MATRIX<T>& matrix)
{
	if (matrix.rows() != matrix.cols())
		throw std::logic_error("Determinant for non-square matrix does not exist.");

	MATRIX<T> tmp = matrix;
	row_echelon_form(tmp);
	return diagonal_elements_product(tmp);
}

template <typename T>
T MATRIX_SoLE<T>::max_norm(const MATRIX<T>& matrix)
{
	T ret = 0;
	for (unsigned i = 0; i < matrix.rows(); i++)
		for (unsigned j = 0; j < matrix.cols(); j++)
			if (ret < abs(matrix.field(i,j)))
				ret = abs(matrix.field(i,j));
	return ret;
}

