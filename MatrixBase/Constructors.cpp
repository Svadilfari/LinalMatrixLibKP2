#include "MatrixClasses.h"

Matrix::Matrix() : N(1), M(1) {
	vector<double> v1;
	v1.push_back(0);
	v.push_back(v1);
}

Matrix::Matrix(const unsigned int& N, const unsigned int& M) : N(N), M(M) {
	if ((N == 0) || (M == 0)) throw ZeroMN();
	v.reserve(M);
	for (unsigned int i = 0; i < M; i++) {
		vector<double> v1;
		v1.reserve(N);
		for (unsigned int j = 0; j < N; j++) v1.push_back(0);
		v.push_back(v1);
	}
};

Matrix::Matrix(const vector<double>& v1) {
	if (v1.size() == 0) throw ZeroMN();
	v.push_back(v1);
	M = 1;
	N = v1.size();
}

Matrix::Matrix(const vector<vector<double>>& v) : v(v) {
	M = v.size();
	if (M == 0) throw ZeroMN();
	N = v[0].size();
	if (N == 0) throw ZeroMN();
	for (auto it = v.begin() + 1; it != v.end(); it++) {
		if (it->size() != N) throw RowsAndCols();
	}
};

vector<vector<double>> IdMatrix::makeId(const unsigned int& num) const {
	if (num == 0) throw ZeroMN();
	vector<vector<double>> rows;
	rows.reserve(num);
	for (int i = 0; i < num; i++) {
		vector<double> v1;
		v1.reserve(num);
		for (int j = 0; j < num; j++) {
			if (i == j) v1.push_back(1);
			else v1.push_back(0);
		}
		rows.push_back(v1);
	}
	return rows;
}

Matrix UpperTriMatrix::makeUpTri(const Matrix& A) {
	if (!A.is_square()) throw NotSquare();
	return forwardGauss(A);
}

vector<vector<double>> UpperTriMatrix::makeFromJagged(vector<vector<double>> v1) {
	if (v1.size() == 0) throw ZeroMN();
	if (v1[0].size() == 0) throw ZeroMN();
	if (v1.size() != v1[0].size()) throw NotSquare();
	vector<vector<double>> v3;
	v3.reserve(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		if (v1[i].size() != v1.size() - i) throw RowsAndCols();
		vector<double> v2;
		v2.reserve(i);
		for (int j = 0; j < i; j++) v2.push_back(0);
		v2.insert(v2.end(), v1[i].begin(), v1[i].end());
		v3.push_back(v2);
	}
	return v3;
}

vector<vector<double>> UpperTriMatrix::parse_utm(const string& fname) {
	vector<vector<double>> v3;
	if (fname.substr(fname.find_last_of(".") + 1) == "txt") {
		ifstream in;
		in.open(fname);
		if (!in.is_open()) throw OpenProblem();
		string line;
		while (getline(in, line)) {
			stringstream stream(line);
			string expr;
			vector<double> v1;
			while (getline(stream, expr, '\t')) {
				if (RUS != 0) if (expr.find(',') != string::npos) expr[expr.find(',')] = '.';
				v1.push_back(atof(expr.c_str()));
			}
			v3.push_back(v1);
		}
		in.close();
	}
	else if (fname.substr(fname.find_last_of(".") + 1) == "bin") {
		ifstream in;
		in.open(fname);
		if (!in.is_open()) throw OpenProblem();
		unsigned int M;
		unsigned int N;
		in.read((char*)&M, sizeof(unsigned int));
		in.read((char*)&N, sizeof(unsigned int));
		if (M != N) throw NotSquare();
		v3.reserve(M);
		for (int i = 0; i < M; i++) {
			vector<double> v1;
			v1.reserve(N);
			for (int j = i; j < N; j++) {
				double elm;
				in.read((char*)&elm, sizeof(double));
				if (in.eof()) throw FormatProblem();
				v1.push_back(elm);
			}
			v3.push_back(v1);
		}
		in.close();
	}
	else {
		throw UnknownExtention();
	}
	return v3;
}

Matrix LowerTriMatrix::makeLowTri(const Matrix& A) {
	if (!A.is_square()) throw NotSquare();
	return backwardGauss(A);
}

vector<vector<double>> LowerTriMatrix::makeFromJagged(vector<vector<double>> v1) {
	if (v1.size() == 0) throw ZeroMN();
	if (v1[0].size() == 0) throw ZeroMN();
	if (v1.size() != v1[v1.size() - 1].size()) throw NotSquare();
	vector<vector<double>> v3;
	v3.reserve(v1.size());
	for (int i = 0; i < v1.size(); i++) {
		if (v1[i].size() != i + 1) throw RowsAndCols();
		vector<double> v2;
		vector<double> v2n5;
		v2n5.reserve(v1.size() - i - 1);
		for (int j = v1.size() - i - 1; j > 0; j--) v2n5.push_back(0);
		v2.insert(v2.end(), v1[i].begin(), v1[i].end());
		v2.insert(v2.end(), v2n5.begin(), v2n5.end());
		v3.push_back(v2);
	}
	return v3;
}

vector<vector<double>> LowerTriMatrix::parse_ltm(const string& fname) {
	vector<vector<double>> v3;
	if (fname.substr(fname.find_last_of(".") + 1) == "txt") {
		ifstream in;
		in.open(fname);
		if (!in.is_open()) throw OpenProblem();
		string line;
		while (getline(in, line)) {
			stringstream stream(line);
			string expr;
			vector<double> v1;
			while (getline(stream, expr, '\t')) {
				if (RUS != 0) if (expr.find(',') != string::npos) expr[expr.find(',')] = '.';
				v1.push_back(atof(expr.c_str()));
			}
			v3.push_back(v1);
		}
		in.close();
	}
	else if (fname.substr(fname.find_last_of(".") + 1) == "bin") {
		ifstream in;
		in.open(fname);
		if (!in.is_open()) throw OpenProblem();
		unsigned int M;
		unsigned int N;
		in.read((char*)&M, sizeof(unsigned int));
		in.read((char*)&N, sizeof(unsigned int));
		if (M != N) throw NotSquare();
		v3.reserve(M);
		for (int i = 0; i < M; i++) {
			vector<double> v1;
			v1.reserve(N);
			for (int j = 0; j < i + 1; j++) {
				double elm;
				in.read((char*)&elm, sizeof(double));
				if (in.eof()) throw FormatProblem();
				v1.push_back(elm);
			}
			v3.push_back(v1);
		}
		in.close();
	}
	else {
		throw UnknownExtention();
	}
	return v3;
}

Matrix DiagonalMatrix::makeDiag(const Matrix& A) {
	if (!A.is_square()) throw NotSquare();
	vector<double> v2;
	v2.reserve(A.get_matrix().size());
	for (int i = 0; i < A.get_matrix().size(); i++) {
		v2.push_back((A.get_matrix())[i][i]);
	}
	return makeFromJagged(v2);
}

vector<vector<double>> DiagonalMatrix::makeFromJagged(vector<double> v1) {
	if (v1.size() == 0) throw ZeroMN();
	Matrix Id1(v1.size(), v1.size());
	vector<vector<double>> v3 = Id1.get_matrix();
	for (int i = 0; i < v1.size(); i++) {
		v3[i][i] = v1[i];
	}
	return v3;
}

vector<double> DiagonalMatrix::parse_diag(const string& fname) const {
	vector<double> v1;
	if (fname.substr(fname.find_last_of(".") + 1) == "txt") {
		ifstream in;
		in.open(fname);
		if (!in.is_open()) throw OpenProblem();
		string expr;
		while (getline(in, expr, '\t')) {
			if (RUS != 0) if (expr.find(',') != string::npos) expr[expr.find(',')] = '.';
			v1.push_back(atof(expr.c_str()));
		}
		in.close();
	}
	else if (fname.substr(fname.find_last_of(".") + 1) == "bin") {
		ifstream in;
		in.open(fname);
		if (!in.is_open()) throw OpenProblem();
		while (!in.eof()) {
			double elm;
			in.read((char*)&elm, sizeof(double));
			v1.push_back(elm);
		}
		in.close();
	}
	else {
		throw UnknownExtention();
	}
	return v1;
}

Matrix SymmetricMatrix::makeSym(vector<vector<double>>& v1) {
	if (v1.size() == 0) throw ZeroMN();
	if (v1.size() == v1[0].size()) {
		UpperTriMatrix A(v1);
		vector<double> v2;
		v2.reserve(v1.size());
		for (auto it = v1.begin(); it != v1.end(); it++) {
			v2.push_back((*it)[0]);
		}
		DiagonalMatrix B(v2);
		return A + A.Transpose() - B;
	}
	else {
		LowerTriMatrix A(v1);
		vector<double> v2;
		v2.reserve(v1.size());
		for (auto it = v1.begin(); it != v1.end(); it++) {
			v2.push_back((*it)[it->size() - 1]);
		}
		DiagonalMatrix B(v2);
		return A + A.Transpose() - B;
	}
};

vector<vector<double>> SymmetricMatrix::parse_sym(const string& fname, bool up) {
	vector<vector<double>> v3;
	if (up) {
		UpperTriMatrix A(fname);
		v3 = (A + A.Transpose() - DiagonalMatrix(A)).get_matrix();
	}
	else {
		LowerTriMatrix A(fname);
		v3 = (A + A.Transpose() - DiagonalMatrix(A)).get_matrix();
	}
	return v3;
}
