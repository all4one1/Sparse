#pragma once
#include <vector>

struct SparseMatrix
{
	enum side
	{
		center,
		east,
		north,
		west,
		south
	};
	// format of storage:
	// Compressed Sparse Row
	int n = 0;

	int nval = 0;
	int ncol = 0;
	int nraw = 0;

	std::vector <double> val;
	std::vector <int> col;
	std::vector <int> raw;
	std::vector <double> diag;
	std::vector <int> type;

	SparseMatrix()
	{
		raw.push_back(0);
	}
	SparseMatrix(int n_)
	{
		resize(n_);
	}
	~SparseMatrix()
	{
		//printf("SparseMatrix destructor executed \n");
	}

	void add(double v, int i, int t = -1)
	{
		type.push_back(t);
		val.push_back(v);
		col.push_back(i);
		nraw++;
	}
	void endline(int l)
	{
		raw[l] = nraw;
	}
	void resize(int n_);

	void make_sparse(int nx, int ny, double a = 1);

	void make_sparse(int n_, double** M);

	void print_index_ij(int l);

	double get_element(int ii, int jj);
	int get_index(int ii, int jj);
	int get_type(int ii, int jj);
	void set_type(int ii, int jj, int t);

	void update_diag();

	double& operator()(int ii, int jj);

	void show_storage();

	void recover_full(int);
	void recover_full2();
	void recover_type();
	void print_all();

	void sequaent_print();

	void erase_zeros();

	// Aij*yj
	double line(int q, double* y);

	double line1(int q, double* y);

	double line2(int q, double* y);

	double max_element_abs();

	double conditionNumber();
	double conditionNumber2();
};
