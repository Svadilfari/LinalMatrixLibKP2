#pragma once
#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<fstream>
#include<sstream>
#include<algorithm>
#include "Problems.h"
#include "FileProblems.h"
using namespace std;
#define ACCURACY 0.0000001
#define RUS 1

#define LAB 2

//Матрица M
class Matrix {
protected:
	vector<vector<double>> v;
	unsigned int M;
	unsigned int N;

	vector<double> addvectors(const vector<double>& v1, const vector<double>& v2) const;
	vector<double> subbvectors(const vector<double>& v1, const vector<double>& v2) const;
	double multiplyvectors(const vector<double>& v1, const vector<double>& v2) const;
	vector<double> amultvector(const double& a, const vector<double>& v1) const;
	double MultDiag(const Matrix& A) const;

public:

	unsigned int get_rows() const { return M; }
	unsigned int get_columns() const { return N; }
	vector<vector<double>> get_matrix() const { return v; }

	bool is_square() const {
		if (M == N) return true;
		else return false;
	};
	bool is_same_size(const Matrix& B) const {
		if ((M == B.M) && (N == B.N)) return true;
		else return false;
	};
	bool is_vector() const {
		if ((N == 1) || (M == 1)) return true;
		else return false;
	};

	Matrix concatenate(const Matrix& B) const;
	double Det() const;
	Matrix Inverse() const;
	Matrix Adamar(const Matrix& B) const;
	double Trace() const;
	double Norm() const;
	double MaxNorm() const;
	int Rank() const;
	Matrix Transpose() const;

	Matrix();
	Matrix(const unsigned int& N, const unsigned int& M);
	Matrix(const vector<double>& v1);
	Matrix(const vector<vector<double>>& v);
	Matrix(const string& fname);
	Matrix(const Matrix& A) : M(A.M), N(A.N), v(A.v){}

	friend Matrix forwardGauss(const Matrix& A);
	friend Matrix backwardGauss(const Matrix& A);
	friend double scalarMult(const Matrix& A, const Matrix& B);
	friend double AngleVectors(const Matrix& A, const Matrix& B);
	friend Matrix operator+(const Matrix& A, const Matrix& B);
	friend Matrix operator-(const Matrix& A, const Matrix& B);
	friend Matrix operator*(const Matrix& A, const Matrix& B);
	friend Matrix operator*(const double& a, const Matrix& A);
	friend Matrix operator*(const Matrix& A, const double& a);
	friend ostream& operator<<(ostream& out, const Matrix& A);
	friend istream& operator>>(istream& in, const Matrix& A);
};

class IdMatrix : public Matrix{
private:
	vector<vector<double>> makeId(const unsigned int& num) const;

public:
	IdMatrix() : Matrix() {};
	IdMatrix(const int& num) : Matrix(makeId(num)) {};
	IdMatrix(const Matrix& A) : Matrix(makeId(A.get_rows())) {};
	IdMatrix(const string& fname) : Matrix(makeId(Matrix(fname).get_rows())){};
};

class UpperTriMatrix : public Matrix {
private:
	Matrix makeUpTri(const Matrix& A);
	vector<vector<double>> makeFromJagged(vector<vector<double>> v1);
	vector<vector<double>> parse_utm(const string& fname);

public:
	UpperTriMatrix() : Matrix() {};
	UpperTriMatrix(const Matrix& A) : Matrix(makeUpTri(A)) {};
	UpperTriMatrix(const vector<vector<double>>& v1) : Matrix(makeFromJagged(v1)) {};
	UpperTriMatrix(const string& fname) : Matrix(makeFromJagged(parse_utm(fname))) {};
};

class LowerTriMatrix : public Matrix {
private:
	Matrix makeLowTri(const Matrix& A);
	vector<vector<double>> makeFromJagged(vector<vector<double>> v1);
	vector<vector<double>> parse_ltm(const string& fname);

public:
	LowerTriMatrix() : Matrix() {};
	LowerTriMatrix(const Matrix& A) : Matrix(makeLowTri(A)) {};
	LowerTriMatrix(const vector<vector<double>>& v1) : Matrix(makeFromJagged(v1)) {};
	LowerTriMatrix(const string& fname) : Matrix(makeFromJagged(parse_ltm(fname))) {};
};

class DiagonalMatrix : public Matrix {
private:
	Matrix makeDiag(const Matrix& A);
	vector<vector<double>> makeFromJagged(vector<double> v1);
	vector<double> parse_diag(const string& fname) const;

public:
	DiagonalMatrix() : Matrix() {};
	DiagonalMatrix(const Matrix& A) : Matrix(makeDiag(A)) {};
	DiagonalMatrix(const vector<double>& v1) : Matrix(makeFromJagged(v1)) {};
	DiagonalMatrix(const string& fname) : Matrix(makeFromJagged(parse_diag(fname))) {};
};

class SymmetricMatrix : public Matrix {
private:
	Matrix makeSym(vector<vector<double>>& v1);
	vector<vector<double>> parse_sym(const string& fname, bool up);

public:
	SymmetricMatrix() : Matrix() {};
	SymmetricMatrix(vector<vector<double>>& v1) : Matrix(makeSym(v1)) {};
	SymmetricMatrix(const string& fname, const bool isUpper = true) : Matrix(parse_sym(fname, isUpper)) {};
};
