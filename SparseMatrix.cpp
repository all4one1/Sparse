#include "SparseMatrix.h"
#include <iostream>
#include <iomanip>
#include <fstream>
using std::cout;
using std::endl;
using std::ofstream;

void SparseMatrix::resize(int n_)
{
	n = n_;
	val.clear(); 	col.clear(); raw.clear();	diag.clear();
	raw.resize(n + 1);	nraw = 0;
	update_diag();
}

void SparseMatrix::make_sparse(int nx, int ny, double a)
{
	int l = -1;
	int off = nx;
	n = nx * ny;
	raw.resize(n + 1);
	raw[0] = 0;
	double ap, ae, aw, an, as;
	ap = 1 + 4 * a;
	ae = -a;
	aw = -a;
	an = -a;
	as = -a;



	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			l = i + off * j;

			ap = 1 + 4 * a;
			if (i == 0) ap -= a;
			if (i == nx - 1) ap -= a;
			if (j == 0) ap -= a;
			if (j == ny - 1) ap -= a;


			if (j > 0)


				add(as, l - off, south);

			if (i > 0)
				add(ae, l - 1, east);

			add(ap, l, center);

			if (i < nx - 1)
				add(aw, l + 1, west);

			if (j < ny - 1)
				add(an, l + off, north);

			endline(l + 1);

		}
	}
	nval = ncol = (int)val.size();

	update_diag();
}

void SparseMatrix::make_sparse(int n_, double** M)
{
	n = n_;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
		{
			if (M[i][j] != 0)
			{
				val.push_back(M[i][j]);
				col.push_back(j);
				nraw++;
			}
			if (i == j)
				diag.push_back(M[i][j]);
		}
		raw.push_back(nraw);
	}
	nval = ncol = (int)val.size();
}

void SparseMatrix::print_index_ij(int l)
{
	if (!(l < nval)) { cout << " out of range " << endl; return; }

	int j = col[l];
	int i;
	for (int q = 0; q < n; q++)
	{
		if (l < raw[q + 1])
		{
			i = q;
			break;
		}
	}

	cout << "[" << i << "][" << j << "]" << endl;
}

double SparseMatrix::get_element(int ii, int jj)
{
	double v = 0;
	for (int j = raw[ii]; j < raw[ii + 1]; j++)
	{
		if (col[j] == jj) v = val[j];
	}
	return v;
}

int SparseMatrix::get_index(int ii, int jj)
{
	int id = -1;
	for (int j = raw[ii]; j < raw[ii + 1]; j++)
	{
		if (col[j] == jj) id = j;
	}
	return id;
}
void SparseMatrix::set_type(int ii, int jj, int t)
{
	int id = -1;
	for (int j = raw[ii]; j < raw[ii + 1]; j++)
	{
		if (col[j] == jj) id = j;
	}
	type[id] = t;
}
int SparseMatrix::get_type(int ii, int jj)
{
	int t = -1;
	for (int j = raw[ii]; j < raw[ii + 1]; j++)
	{
		if (col[j] == jj) t = type[j];
	}
	return t;
}
void SparseMatrix::update_diag()
{
	diag.resize(n);
	for (int i = 0; i < n; i++)
	{
		diag[i] = get_element(i, i);
	}
}

double& SparseMatrix::operator()(int ii, int jj)
{
	//cout << "index: " << ii << " " << jj << endl;
	if (!(ii < n) || !(jj < n))
	{
		cout << "bad case " << endl;
		double* v = new double; 	return *v;
	}


	int l1 = raw[ii];
	int l2 = raw[ii + 1];


	if (l2 - l1 == 0)
	{
		val.insert(val.begin() + l1, 0.0);
		col.insert(col.begin() + l1, jj);
		type.insert(type.begin() + l1, -1);
		for (int i = ii + 1; i < (int)raw.size(); i++)
		{
			raw[i]++;
		}
		nval++;
		return val[l1];
	}


	//before 
	if (jj >= 0 && jj < col[l1])
	{
		val.insert(val.begin() + l1, 0.0);
		col.insert(col.begin() + l1, jj);
		type.insert(type.begin() + l1, -1);
		for (int i = ii + 1; i < (int)raw.size(); i++)
		{
			raw[i]++;
		}
		nval++;
		return val[l1];
	}

	//after
	if (jj < n && jj > col[l2 - 1])
	{
		val.insert(val.begin() + l2, 0.0);
		col.insert(col.begin() + l2, jj);
		type.insert(type.begin() + l2, -1);
		for (int i = ii + 1; i < (int)raw.size(); i++)
		{
			raw[i]++;
		}
		nval++;
		return val[l2];
	}
	//inside existing
	for (int j = raw[ii]; j < raw[ii + 1]; j++)
	{
		if (col[j] == jj)
		{
			return val[j];
		}
	}


	//inside non-existing
	for (int j = raw[ii]; j < raw[ii + 1]; j++)
	{
		if (jj > col[j] && jj < col[j + 1])
		{
			val.insert(val.begin() + j + 1, 0.0);
			col.insert(col.begin() + j + 1, jj);
			type.insert(type.begin() + j + 1, -1);
			for (int i = ii + 1; i < (int)raw.size(); i++)
			{
				raw[i]++;
			}
			nval++;
			return val[j + 1];

		}

	}

	/*
	//inside non-existing
	for (int j = raw[ii]; j < raw[ii + 1]; j++)
	{
		for (int q = col[j] + 1; q < col[j + 1]; q++)
		{
			if (jj == q)
			{
				//cout << "jj:" << j << " " << q << endl;
				val.insert(val.begin() + j, 0.0);
				col.insert(col.begin() + j, jj);
				//cout << "q:" << q << endl;

				for (int i = ii + 1; i < (int)raw.size(); i++)
				{
					raw[i]++;
				}
				nval++;
				return val[j];
			}
		}

	}
	*/

	cout << "bad case-end " << endl;
	double* v = new double; 	return *v;
}

void SparseMatrix::show_storage()
{
	double S = double(val.capacity() * 8 + col.capacity() * 4 + raw.capacity() * 4 + diag.capacity() * 8 + type.capacity() * 4);
	int num = int(val.capacity());
	cout << num << " elements, " << S / 1024 / 1024 << " MB approx. matrix memory usage" << " \n\n";
}

void SparseMatrix::recover_full(int precision = 4)
{
	ofstream out("coef.dat");
	ofstream out2("info.dat");
	out << std::fixed << std::setprecision(precision);
	for (unsigned int k = 0; k < raw.size() - 1; k++)
	{
		std::vector <double> line(n);
		std::vector <int> t(n);
		for (int i = 0; i < n; i++)
		{
			line[i] = 0.0;
			//	t[i] = -1;
		}

		for (int j = raw[k]; j < raw[k + 1]; j++)
		{
			line[col[j]] = val[j];
			//	t[col[j]] = type[j];
		}

		for (int i = 0; i < n; i++)
		{
			out << line[i] << " ";
			//	if (t[i] == -1) out2 << "_" << " ";
			//	else out2 << t[i] << " ";
		}
		out << endl;
		out2 << endl;
	}
}

void SparseMatrix::recover_full2()
{
	ofstream m("matrix.dat");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
		{
			m << get_element(i, j) << " ";
		}
		m << endl;
	}
}

void SparseMatrix::recover_type()
{
	ofstream m("type.dat");
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
		{
			std::string str = "_____";
			int t = get_type(i, j);
			if (t == center) str = "centr";
			if (t == south) str = "south";
			if (t == north) str = "north";
			if (t == west) str = "west_";
			if (t == east) str = "east_";

			m << str << " ";
		}
		m << endl;
	}
}

void SparseMatrix::print_all()
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
		{
			std::cout << get_element(i, j) << " ";
		}
		std::cout << endl;
	}
}

void SparseMatrix::erase_zeros()
{
	//очистка нулей после LU-разложения нужна тут
}

double SparseMatrix::line(int q, double* y)
{
	double s = 0.0;

	for (int j = raw[q]; j < raw[q + 1]; j++)
	{
		s += val[j] * y[col[j]];
	}
	return s;
}

double SparseMatrix::line1(int q, double* y)
{
	double s = 0.0;
	for (int j = raw[q]; j < raw[q + 1]; j++)
	{
		if (col[j] >= q) break;
		s += val[j] * y[col[j]];
	}
	return s;
}

double SparseMatrix::line2(int q, double* y)
{
	double s = 0.0;
	for (int j = raw[q]; j < raw[q + 1]; j++)
		//for (int j = raw[q + 1] - 1; j <= raw[q]; j--)
	{
		if (col[j] >= q + 1) {
			s += val[j] * y[col[j]];
		}
	}
	return s;
}

void SparseMatrix::sequaent_print()
{
	ofstream seq("seq.dat");
	seq << nval << endl;
	for (int i = 0; i < nval; i++)
		seq << val[i] << " ";
	seq << endl;
	for (int i = 0; i < nval; i++)
		seq << col[i] << " ";
	seq << endl;

	for (int i = 0; i < n + 1; i++)
		seq << raw[i] << " ";
}

double SparseMatrix::max_element_abs()
{
	double max = 0;
	double v = 0;
	for (int i = 0; i < nval; i++)
	{
		v = abs(val[i]);
		if (v > max)
			max = v;
	}

	return max;
}