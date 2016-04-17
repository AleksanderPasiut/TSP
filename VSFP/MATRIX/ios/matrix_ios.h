// This file provides additional functions for MATRIX template class to be
// comfortably used with standard input/output library (iostream).

#pragma once

#include <iostream>
#include <iomanip>
using namespace std;

#include "..\matrix.h"

template <typename T>
istream& operator>> (istream& str, MATRIX<T>& matrix)
{
	for (unsigned i = 0; i < matrix.rows(); i++)
		for (unsigned j = 0; j < matrix.cols(); j++)
			str >> matrix.field(i,j);

	return str;
}

template <typename T>
ostream& operator<< (ostream& str, MATRIX<T>& matrix)
{
	for (unsigned i = 0; i < matrix.rows(); i++)
	{
		for (unsigned j = 0; j < matrix.cols(); j++)
			str << setw(6) << matrix.field(i,j) << " ";

		str << endl;
	}

	return str;
}
