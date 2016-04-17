#pragma once

#include <stdexcept>

template<typename T>
MATRIX<T>::MATRIX(unsigned rows, unsigned cols)
{
	if (!rows)
		throw std::logic_error("Rows amount is zero.");

	if (!cols)
		throw std::logic_error("Columns amount is zero.");
	
	rows_amount = rows;
	cols_amount = cols;

	table = new T* [rows];
	for (unsigned i = 0; i < rows; i++)
	{
		table[i] = new T[cols];

		for (unsigned j = 0; j < cols; j++)
			table[i][j] = 0;
	}
}

template<typename T>
MATRIX<T>::MATRIX(const MATRIX<T>& arg)
{
	table = new T* [arg.rows_amount];
	for (unsigned i = 0; i < arg.rows_amount; i++)
	{
		table[i] = new T[arg.cols_amount];

		for (unsigned j = 0; j < arg.cols_amount; j++)
			table[i][j] = arg.table[i][j];
	}

	rows_amount = arg.rows_amount;
	cols_amount = arg.cols_amount;
}

template<typename T>
MATRIX<T>::~MATRIX()
{
	for (unsigned i = 0; i < rows_amount; i++)
		delete[] table[i];
	delete[] table;
}

template<typename T>
MATRIX<T>& MATRIX<T>::operator= (const MATRIX<T>& arg)
{
	if (rows_amount == arg.rows_amount && cols_amount == arg.cols_amount)
	{
		for (unsigned i = 0; i < rows_amount; i++)
			for (unsigned j = 0; j < cols_amount; j++)
				table[i][j] = arg.table[i][j];
	}
	else
	{
		T** new_table = new T* [arg.rows_amount];
		for (unsigned i = 0; i < arg.rows_amount; i++)
		{
			new_table[i] = new T [arg.cols_amount];
			for (unsigned j = 0; j < arg.cols_amount; j++)
				new_table[i][j] = arg.table[i][j];
		}

		for (unsigned i = 0; i < rows_amount; i++)
			delete[] table[i];
		delete[] table;

		table = new_table;
		rows_amount = arg.rows_amount;
		cols_amount = arg.cols_amount;
	}

	return *this;
}

template<typename T>
bool MATRIX<T>::operator== (const MATRIX<T>& arg) const
{
	return !(*this != arg);
}

template<typename T>
bool MATRIX<T>::operator!= (const MATRIX<T>& arg) const
{
	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < cols_amount; j++)
			if (table[i][j] != arg.table[i][j])
				return true;

	return false;
}

template<typename T>
MATRIX<T>& MATRIX<T>::operator+= (const MATRIX<T>& arg)
{
	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < cols_amount; j++)
			table[i][j] += arg.table[i][j];

	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::operator-= (const MATRIX<T>& arg)
{
	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < cols_amount; j++)
			table[i][j] -= arg.table[i][j];

	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::operator*= (const MATRIX<T>& arg)
{
	*this = *this * arg;
	return *this;
}

template<typename T>
MATRIX<T> MATRIX<T>::operator+ (const MATRIX<T>& arg) const
{
	MATRIX<T> ret = *this;
	return ret += arg;
}

template<typename T>
MATRIX<T> MATRIX<T>::operator- (const MATRIX<T>& arg) const
{
	MATRIX<T> ret = *this;
	return ret -= arg;
}

template<typename T>
MATRIX<T> MATRIX<T>::operator* (const MATRIX<T>& arg) const
{
	if (cols_amount != arg.rows_amount)
		throw std::logic_error("Matrices sizes don't match.");

	MATRIX<T> ret(rows_amount, arg.cols_amount);

	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < arg.cols_amount; j++)
		{
			ret.table[i][j] = 0;
			for (unsigned k = 0; k < cols_amount; k++)
				ret.table[i][j] += table[i][k]*arg.table[k][j];
		}

	return ret;
}

template<typename T>
MATRIX<T>& MATRIX<T>::operator+= (const T& arg)
{
	for (unsigned i = 0; i < min(rows_amount, cols_amount); i++)
		table[i][i] += arg;

	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::operator-= (const T& arg)
{
	for (unsigned i = 0; i < min(rows_amount, cols_amount); i++)
		table[i][i] -= arg;

	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::operator*= (const T& arg)
{
	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < cols_amount; j++)
			table[i][j] *= arg;

	return *this;
}

template<typename T>
MATRIX<T> MATRIX<T>::operator+ (const T& arg) const
{
	MATRIX<T> ret = *this;
	return ret += arg;
}

template<typename T>
MATRIX<T> MATRIX<T>::operator- (const T& arg) const
{
	MATRIX<T> ret = *this;
	return ret -= arg;
}

template<typename T>
MATRIX<T> MATRIX<T>::operator* (const T& arg) const
{
	MATRIX<T> ret = *this;
	return ret *= arg;
}

template<typename T>
T& MATRIX<T>::field(unsigned row, unsigned col)
{
	if (row >= rows_amount)
		throw std::logic_error("Row index is out of range.");

	if (col >= cols_amount)
		throw std::logic_error("Column index is out of range.");

	return table[row][col];
}

template<typename T>
const T& MATRIX<T>::field(unsigned row, unsigned col) const
{
	if (row >= rows_amount)
		throw std::logic_error("Row index is out of range.");

	if (col >= cols_amount)
		throw std::logic_error("Column index is out of range.");

	return table[row][col];
}

template<typename T>
MATRIX<T>& MATRIX<T>::change_size(unsigned rows, unsigned cols)
{
	if (!rows)
		throw std::logic_error("Rows amount is zero.");

	if (!cols)
		throw std::logic_error("Columns amount is zero.");

	T** new_table = new T* [rows];
	for (unsigned i = 0; i < rows; i++)
	{
		new_table[i] = new T [cols];

		for (unsigned j = 0; j < cols; j++)
			new_table[i][j] = (i < rows_amount && j < cols_amount) ? table[i][j] : 0;
	}

	for (unsigned i = 0; i < rows_amount; i++)
		delete[] table[i];
	delete table;

	table = new_table;
	rows_amount = rows;
	cols_amount = cols;
	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::remove_rows(unsigned first, unsigned amount)
{
	if (!amount)
		throw std::logic_error("Amount is zero.");

	if (!first && amount == rows_amount)
		throw std::logic_error("Amount is equal rows amount.");

	if (first+amount > rows_amount)
		throw std::logic_error("There are not enough rows to delete.");

	T** new_table = new T* [rows_amount-amount];
	unsigned i = 0;
	for (; i < first; i++)
		new_table[i] = table[i];
	for (; i < first+amount; i++)
		delete[] table[i];
	for (; i < rows_amount; i++)
		new_table[i-amount] = table[i];

	delete[] table;
	table = new_table;
	rows_amount -= amount;
	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::remove_cols(unsigned first, unsigned amount)
{
	if (!amount)
		throw std::logic_error("Amount is zero.");

	if (!first && amount == cols_amount)
		throw std::logic_error("Amount is equal cols amount.");

	if (first+amount > cols_amount)
		throw std::logic_error("There are not enough cols to delete.");

	for (unsigned i = 0; i < rows_amount; i++)
	{
		T* new_row = new T[cols_amount-amount];
		unsigned j = 0;
		for (; j < first; j++)
			new_row[j] = table[i][j];
		j += amount;
		for (; j < cols_amount; j++)
			new_row[j-amount] = table[i][j];

		delete[] table[i];
		table[i] = new_row;
	}

	cols_amount -= amount;
	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::row_combine(const MATRIX<T>& arg)
{
	if (cols_amount != arg.cols_amount)
		throw std::logic_error("Matrices have different amounts of columns.");

	unsigned shift = rows_amount;
	change_size(rows_amount+arg.rows_amount, cols_amount);

	for (unsigned i = 0; i < arg.rows_amount; i++)
		for (unsigned j = 0; j < cols_amount; j++)
			table[shift+i][j] = arg.table[i][j];

	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::col_combine(const MATRIX<T>& arg)
{
	if (rows_amount != arg.rows_amount)
		throw std::logic_error("Matrices have different amounts of rows.");

	unsigned shift = cols_amount;
	change_size(rows_amount, cols_amount+arg.cols_amount);

	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < arg.cols_amount; j++)
			table[i][shift+j] = arg.table[i][j];

	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::zero()
{
	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < rows_amount; j++)
			table[i][j] = 0;
	return *this;
}

template<typename T>
MATRIX<T>& MATRIX<T>::identity()
{
	if (rows_amount != cols_amount)
		throw std::logic_error("Rectangular matrix cannot be identity matrix.");

	for (unsigned i = 0; i < rows_amount; i++)
		for (unsigned j = 0; j < cols_amount; j++)
			table[i][j] = (i == j) ? 1 : 0;

	return *this;
}
