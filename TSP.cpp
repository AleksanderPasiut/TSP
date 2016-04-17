#include <iostream>
#include "VSFP\\vsfp.h"
using namespace std;

double minLB = HUGE_VAL;

struct NODE
{
	MATRIX<double> m;
	double LB;
	double beg;
	double end;
	NODE* parent;
	NODE* left;
	NODE* right;

	NODE() : left(0), right(0), parent(0) {}
	~NODE()	{ delete left; delete right; }
};

void FillDefaultMatrix(MATRIX<double>& m)
{
	for (unsigned i = 1; i < m.rows(); i++)
	{
		m.field(0, i) = i;
		m.field(i, 0) = i;
	}

	m.field(1, 2) = 3;
	m.field(1, 3) = 93;
	m.field(1, 4) = 13;
	m.field(1, 5) = 33;
	m.field(1, 6) = 9;

	m.field(2, 1) = 4;
	m.field(2, 3) = 77;
	m.field(2, 4) = 42;
	m.field(2, 5) = 21;
	m.field(2, 6) = 16;

	m.field(3, 1) = 45;
	m.field(3, 2) = 17;
	m.field(3, 4) = 36;
	m.field(3, 5) = 16;
	m.field(3, 6) = 28;

	m.field(4, 1) = 39;
	m.field(4, 2) = 90;
	m.field(4, 3) = 80;
	m.field(4, 5) = 56;
	m.field(4, 6) = 7;

	m.field(5, 1) = 28;
	m.field(5, 2) = 46;
	m.field(5, 3) = 88;
	m.field(5, 4) = 33;
	m.field(5, 6) = 25;

	m.field(6, 1) = 3;
	m.field(6, 2) = 88;
	m.field(6, 3) = 18;
	m.field(6, 4) = 46;
	m.field(6, 5) = 92;
}
void AddTopRowNCol(MATRIX<double>& m)
{
	MATRIX<double> topRow(1, m.cols());
	for (unsigned i = 0; i < m.cols(); i++)
		topRow.field(0, i) = static_cast<double>(i)+1.0;

	m = topRow.row_combine(m);

	MATRIX<double> leftCol(m.rows());
	for (unsigned i = 1; i < m.rows(); i++)
		leftCol.field(i, 0) = static_cast<double>(i);

	m = leftCol.col_combine(m);
}
void SetInfinities(MATRIX<double>& m)
{
	for (unsigned i = 0; i < m.rows(); i++)
		m.field(i,i) = HUGE_VAL;
}

double ReduceMatrix(MATRIX<double>& m)
{
	double LB = 0;

	for (unsigned i = 1; i < m.rows(); i++)
	{
		double factor = HUGE_VAL;

		for (unsigned j = 1; j < m.cols(); j++)
			if (factor > m.field(i,j))
				factor = m.field(i,j);

		if (factor == HUGE_VAL)
			factor = 0;

		LB += factor;

		for (unsigned j = 1; j < m.cols(); j++)
			m.field(i,j) -= factor;
	}

	for (unsigned i = 1; i < m.cols(); i++)
	{
		double factor = HUGE_VAL;

		for (unsigned j = 1; j < m.rows(); j++)
			if (factor > m.field(j,i))
				factor = m.field(j,i);

		if (factor == HUGE_VAL)
			factor = 0;

		LB += factor;

		for (unsigned j = 1; j < m.rows(); j++)
			m.field(j,i) -= factor;
	}

	return LB;
}
bool IsMatrixValid(MATRIX<double>& m)
{
	if (m.rows() <= 3 || m.cols() <= 3)
		return false;

	for (unsigned i = 1; i < m.rows(); i++)
		for (unsigned j = 1; j < m.cols(); j++)
			if (m.field(i,j) != 0 && m.field(i,j) != HUGE_VAL)
				return true;

	return false;
}
double FindDivisionPoint(MATRIX<double>& m, double& beg, double& end)
{
	double ret = 0;
	for (unsigned i = 1; i < m.rows(); i++)
		for (unsigned j = 1; j < m.cols(); j++)
			if (m.field(i,j) == 0)
			{
				double row_val = HUGE_VAL;
				for (unsigned k = 1; k < m.cols(); k++)
					if (k != j && row_val > m.field(i,k))
						row_val = m.field(i,k);

				double col_val = HUGE_VAL;
				for (unsigned k = 1; k < m.rows(); k++)
					if (k != i && col_val > m.field(k,j))
						col_val = m.field(k,j);

				if (row_val != HUGE_VAL && 
					col_val != HUGE_VAL && 
					ret < row_val+col_val)
				{
					ret = row_val+col_val;
					beg = m.field(i, 0);
					end = m.field(0, j);
				}
			}

	return ret;
}
void SetDivisionPointInfinity(NODE& node)
{
	for (unsigned i = 0; i < node.m.rows(); i++)
		if (node.beg == node.m.field(i, 0))
			for (unsigned j = 0; j < node.m.cols(); j++)
				if (node.end == node.m.field(0, j))
				{
					node.right->m.field(i,j) = HUGE_VAL;
					break;
				}
}
void SetInvDivisionPointInfinity(NODE& node)
{
	for (unsigned i = 0; i < node.m.rows(); i++)
		if (node.end == node.m.field(i, 0))
			for (unsigned j = 0; j < node.m.cols(); j++)
				if (node.beg == node.m.field(0, j))
				{
					node.left->m.field(i,j) = HUGE_VAL;
					break;
				}
}
void RemoveDivisionPointRowNCol(NODE& node)
{
	for (unsigned i = 0; i < node.m.rows(); i++)
		if (node.beg == node.m.field(i, 0))
			for (unsigned j = 0; j < node.m.cols(); j++)
				if (node.end == node.m.field(0, j))
				{
					node.left->m.remove_rows(i);
					node.left->m.remove_cols(j);
					break;
				}
}
void PrepareLeftNodeMatrix(NODE& node)
{
	node.left = new NODE;
	node.left->m = node.m;
	node.left->LB = node.LB;
	node.left->parent = &node;
	
	SetInvDivisionPointInfinity(node);
	RemoveDivisionPointRowNCol(node);
}
void PrepareRightNodeMatrix(NODE& node)
{
	node.right = new NODE;
	node.right->m = node.m;
	node.right->LB = node.LB;
	node.right->parent = &node;

	SetDivisionPointInfinity(node);
}

void PerformDivision(NODE& node)
{
	node.LB += ReduceMatrix(node.m);

	if (!IsMatrixValid(node.m) || (node.parent && node.m == node.parent->m))
	{
		if (minLB > node.LB)
			minLB = node.LB;
			
		return;
	}

	if (node.LB > minLB)
	{
		if (node.parent->left == &node)
			node.parent->left = 0;
		else if (node.parent->right == &node)
			node.parent->right = 0;

		delete &node;
		return;
	}

	FindDivisionPoint(node.m, node.beg, node.end);
	PrepareLeftNodeMatrix(node);
	PrepareRightNodeMatrix(node);

	PerformDivision(*node.left);
	PerformDivision(*node.right);
}

void MoveBack(NODE* node)
{
	if (!node->parent)
	{
		cout << endl;
		return;
	}

	MoveBack(node->parent);

	if (node->parent->left == node)
		cout << node->parent->beg << " - " << node->parent->end << endl;
}
void WritePath(NODE* node)
{
	if (node->left)
		WritePath(node->left);
	if (node->right)
		WritePath(node->right);
	if (!node->left && !node->right)
		MoveBack(node);
}
void RetrievePaths(NODE* node)
{
	cout << "LB: " << minLB << endl << endl;

	WritePath(node);
}

int main()
{
	NODE head;

	unsigned size = 6;
	cout << "Enter number of cities: ";
	cin >> size;

	head.m.change_size(size, size);
	cout << "Fill factors' matrix (diagonal elements will be ignored): " << endl;
	cin >> head.m;

	//FillDefaultMatrix(head.m);
	AddTopRowNCol(head.m);
	SetInfinities(head.m);
	head.LB = 0;

	cout << head.m << endl;

	PerformDivision(head);

	RetrievePaths(&head);

	cin.get();
	cin.ignore();
	return 0;
}