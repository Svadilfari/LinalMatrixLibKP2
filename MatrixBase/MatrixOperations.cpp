#include "Matrixclasses.h"

vector<double> Matrix::addvectors(const vector<double>& v1, const vector<double>& v2) const {
	vector<double> v3;
	auto it2 = v2.begin();
	v3.reserve(v1.size());
	for (auto it = v1.begin(); it != v1.end(); it++, it2++) {
		v3.push_back(*it + *it2);
	}
	return v3;
};

vector<double> Matrix::subbvectors(const vector<double>& v1, const vector<double>& v2) const {
	vector<double> v3;
	auto it2 = v2.begin();
	v3.reserve(v1.size());
	for (auto it = v1.begin(); it != v1.end(); it++, it2++) {
		v3.push_back(*it - *it2);
	}
	return v3;
};

double Matrix::multiplyvectors(const vector<double>& v1, const vector<double>& v2) const {
	double scal = 0;
	auto it2 = v2.begin();
	for (auto it = v1.begin(); it != v1.end(); it++, it2++) {
		scal += *it * *it2;
	}
	return (scal);
};

vector<double> Matrix::amultvector(const double& a, const vector<double>& v1) const {
	vector<double> v3;
	v3.reserve(v1.size());
	for (auto it = v1.begin(); it != v1.end(); it++) v3.push_back(*it * a);
	return v3;
};

double Matrix::MultDiag(const Matrix& A) const {
	double det = 1;
	for (int i = 0; i < A.N; i++) det *= A.v[i][i];
	return det;
}

double Matrix::Det() const {
	if (!is_square()) throw NotSquare();
	vector<vector<double>> rows = this->v;
	double mv_vec = 1;
	for (int j = 0; j != rows.size(); j++) {
		int st = -1;
		int i = 0;
		for (auto it = rows.begin() + j; it != rows.end(); it++) {
			if (st < 0) {
				if (fabs((*it)[j]) < ACCURACY) {
					i++;
					continue;
				}
				st = i + j;
				if (i > 0) {
					vector<double> tmp = rows[j];
					rows[j] = rows[st];
					rows[st] = tmp;
					mv_vec *= -1;
				}
				continue;
			}
			if (fabs((*it)[j]) < ACCURACY) continue;
			double koef = (*it)[j] / rows[j][j];
			vector<double> tosubb = std::move(this->amultvector(koef, rows[j]));
			vector<double> subbed = std::move(this->subbvectors(*it, tosubb));
			*it = std::move(subbed);
		}
		if (st < 0) rows = std::move(Matrix(this->M, this->N).get_matrix());
	}
	return mv_vec * MultDiag(Matrix(rows));
}

Matrix Matrix::Inverse() const {
	if (!is_square()) throw NotSquare();
	if (fabs(Det()) < ACCURACY) throw NoInverse();
	IdMatrix B(M);
	Matrix C = backwardGauss(forwardGauss(this->concatenate(B)));
	vector<vector<double>> v3;
	v3.reserve(C.v.size());
	int i = 0;
	for (auto it = C.v.begin(); it != C.v.end(); it++, i++) {
		vector<double> v1;
		v1.insert(v1.end(), it->begin() + M, it->end());
		double koef = 1 / (*it)[i];
		v3.push_back(amultvector(koef, v1));
	}
	return Matrix(v3);
}

Matrix Matrix::concatenate(const Matrix& B) const {
	if (M != B.M) throw ConcatProblem();
	vector<vector<double>> v3 = v;
	auto it2 = B.v.begin();
	for (auto it = v3.begin(); it != v3.end(); it++) {
		it->insert(it->end(), it2->begin(), it2->end());
		it2++;
	}
	return Matrix(v3);
}

Matrix Matrix::Adamar(const Matrix& B) const{
	if (!is_same_size(B)) throw SizeProblem();
	vector<vector<double>> v3;
	v3.reserve(v.size());
	auto it2 = B.v.begin();
	for (auto it = v.begin(); it != v.end(); it++, it2++) {
		vector<double> v1;
		v1.reserve(it->size());
		auto j = it2->begin();
		for (auto i = it->begin(); i != it->end(); i++, j++) {
			v1.push_back(*i * *j);
		}
		v3.push_back(v1);
	}
	return Matrix(v3);
}

double Matrix::Trace() const {
	int a = (M < N) ? M : N;
	double tr = 0;
	for (int i = 0; i < a; i++) {
		tr += v[i][i];
	}
	return tr;
}

double Matrix::Norm() const {
	if (is_vector()) {
		return pow(scalarMult(*this, *this), 0.5);
	}
	else {
		double norm = 0;
		Matrix A = Adamar(*this);
		for (auto it = A.v.begin(); it != A.v.end(); it++) {
			for (auto j = it->begin(); j != it->end(); j++) {
				norm += *j;
			}
		}
		return pow(norm, 0.5);
	}
}

double Matrix::MaxNorm() const{
	if (!is_vector()) throw NotAVector();
	vector<vector<double>> v3;
	if (N == 1)  v3 = this->Transpose().get_matrix();
	else v3 = get_matrix();
	return *std::max_element(v3[0].begin(), v3[0].end());
}

double AngleVectors(const Matrix& A, const Matrix& B) {
	if (!A.is_vector()) throw NotAVector();
	if (!B.is_vector()) throw NotAVector();
	if ((A.M * A.N) != (B.M * B.N)) throw SizeProblem();
	double a = A.Norm();
	double b = B.Norm();
	if ((a < ACCURACY) || (b < ACCURACY)) return 0;
	double c = scalarMult(A, B);
	return acos(c / (a * b));
}

Matrix forwardGauss(const Matrix& A) {
	vector<vector<double>> rows = A.v;
	int Nmin = (A.M < A.N) ? A.M : A.N;
	for (int j = 0; j != Nmin; j++) {
		int st = -1;
		int i = 0;
		for (auto it = rows.begin() + j; it != rows.end(); it++) {
			if (st < 0) {
				if (fabs((*it)[j]) < ACCURACY) {
					i++;
					continue;
				}
				st = i + j;
				if (i > 0) {
					vector<double> tmp = rows[j];
					rows[j] = rows[st];
					rows[st] = tmp;
				}
				continue;
			}
			if (fabs((*it)[j]) < ACCURACY) continue;
			double koef = (*it)[j] / rows[j][j];
			vector<double> tosubb = std::move(A.amultvector(koef, rows[j]));
			vector<double> subbed = std::move(A.subbvectors(*it, tosubb));
			*it = std::move(subbed);
		}
	}
	return Matrix(rows);
}

Matrix backwardGauss(const Matrix& A) {
	vector<vector<double>> rows = A.v;
	int N = rows.size() - 1;
	for (int j = 0; j != rows.size(); j++) {
		int st = -1;
		int i = 0;
		for (auto it = rows.rbegin() + j; it != rows.rend(); it++) {
			if (st < 0) {
				if (fabs((*it)[N - j]) < ACCURACY) {
					i++;
					continue;
				}
				st = N - i - j;
				if (i > 0) {
					vector<double> tmp = rows[N - j];
					rows[N - j] = rows[st];
					rows[st] = tmp;
				}
				continue;
			}
			if (fabs((*it)[N - j]) < ACCURACY) continue;
			double koef = (*it)[N - j] / rows[N - j][N - j];
			vector<double> tosubb = std::move(A.amultvector(koef, rows[N - j]));
			vector<double> subbed = std::move(A.subbvectors(*it, tosubb));
			*it = std::move(subbed);
		}
	}
	return Matrix(rows);
}

int Matrix::Rank() const {
	vector<vector<double>> rows;
	if (N > M) rows = (this->Transpose()).get_matrix();
	else rows = this->v;
	int Nmin = (this->M < this->N) ? M : N;
	for (int j = 0; j != Nmin; j++) {
		int st = -1;
		int i = 0;
		for (auto it = rows.begin() + j; it != rows.end(); it++) {
			if (st < 0) {
				if (fabs((*it)[j]) < ACCURACY) {
					i++;
					continue;
				}
				st = i + j;
				if (i > 0) {
					vector<double> tmp = rows[j];
					rows[j] = rows[st];
					rows[st] = tmp;
				}
				continue;
			}
			if (fabs((*it)[j]) < ACCURACY) continue;
			double koef = (*it)[j] / rows[j][j];
			vector<double> tosubb = std::move(this->amultvector(koef, rows[j]));
			vector<double> subbed = std::move(this->subbvectors(*it, tosubb));
			*it = std::move(subbed);
		}
	}
	for (int i = Nmin - 1; i >= 0; i--) {
		int st2 = -1;
		for (int j = i; j >= 0; j--) {
			if (st2 == -1) {
				if (fabs(rows[j][i]) < ACCURACY) continue;
				else st2 = j;
				vector<double> tmp = rows[i];
				rows[i] = rows[st2];
				rows[st2] = tmp;
				continue;
			}
			if (fabs(rows[j][i]) > ACCURACY) {
				double koef = rows[j][i] / rows[i][i];
				vector<double> tosubb = std::move(this->amultvector(koef, rows[i]));
				vector<double> subbed = std::move(this->subbvectors(rows[j], tosubb));
				rows[j] = std::move(subbed);
			}
		}
	}
	int r = 0;
	for (auto it = rows.begin(); it != rows.begin() + Nmin; it++) {
		for (auto j = it->begin(); j != it->end(); j++) {
			if (fabs(*j) < ACCURACY) continue;
			else {
				r++;
				break;
			}
		}
	}
	return r;
}

Matrix Matrix::Transpose() const {
	vector<vector<double>> v3;
	v3.reserve(v[0].size());
	int i = 0;
	for (auto it = v[0].begin(); it != v[0].end(); it++, i++) {
		vector<double> v1;
		v1.reserve(v.size());
		for (auto it2 = v.begin(); it2 != v.end(); it2++) {
			v1.push_back((*it2)[i]);
		}
		v3.push_back(v1);
	}
	return Matrix(v3);
}

double scalarMult(const Matrix& A, const Matrix& B) {
	if (!A.is_vector()) throw NotAVector();
	if (!B.is_vector()) throw NotAVector();
	if ((A.M * A.N) != (B.M * B.N)) throw SizeProblem();
	vector<double> v1;
	vector<double> v2;
	if (A.M != 1) {
		v1.reserve(A.v.size());
		for (auto it = A.v.begin(); it != A.v.end(); it++) {
			v1.push_back((*it)[0]);
		}
	}
	else v1 = A.v[0];
	if (B.M != 1) {
		v2.reserve(B.v.size());
		for (auto it = B.v.begin(); it != B.v.end(); it++) {
			v2.push_back((*it)[0]);
		}
	}
	else v2 = B.v[0];
	return A.multiplyvectors(v1, v2);
}

Matrix operator+(const Matrix& A, const Matrix& B) {
	if (!A.is_same_size(B)) throw SizeProblem();
	vector<vector<double>> v3;
	v3.reserve(A.v.size());
	auto it2 = B.v.begin();
	for (auto it = A.v.begin(); it != A.v.end(); it++, it2++) {
		v3.push_back(A.addvectors(*it, *it2));
	}
	return Matrix(v3);
}

Matrix operator-(const Matrix& A, const Matrix& B) {
	if (!A.is_same_size(B)) throw SizeProblem();
	vector<vector<double>> v3;
	v3.reserve(A.v.size());
	auto it2 = B.v.begin();
	for (auto it = A.v.begin(); it != A.v.end(); it++, it2++) {
		v3.push_back(A.subbvectors(*it, *it2));
	}
	return Matrix(v3);
}

Matrix operator*(const Matrix& A, const Matrix& B) {
	if (!A.is_same_size(B)) throw SizeProblem();
	vector<vector<double>> v3;
	v3.reserve(A.v.size());
	for (auto it = A.v.begin(); it != A.v.end(); it++) {
		vector<double> v1;
		v1.reserve(B.N);
		for (int j = 0; j < B.N; j++) {
			vector<double> v2;
			v2.reserve(B.v.size());
			for (auto it2 = B.v.begin(); it2 != B.v.end(); it2++) {
				v2.push_back((*it2)[j]);
			}
			v1.push_back(A.multiplyvectors(*it, v2));
		}
		v3.push_back(v1);
	}
	return Matrix(v3);
}

Matrix operator*(const double& a, const Matrix& A) {
	vector<vector<double>> v3;
	v3.reserve(A.v.size());
	for (auto it = A.v.begin(); it != A.v.end(); it++) {
		vector<double> v1;
		v1.reserve(it->size());
		for (auto j = it->begin(); j != it->end(); j++) {
			v1.push_back(*j * a);
		}
		v3.push_back(v1);
	}
	return Matrix(v3);
}

Matrix operator*(const Matrix& A, const double& a) {
	return a * A;
}

ostream& operator<< (ostream& out, const Matrix& A){
	if (RUS != 0) out.imbue(locale("rus_rus.1251"));
	for (auto it = A.v.begin(); it != A.v.end(); it++) {
		for (auto j = it->begin(); j != it->end(); j++) {
			out << ((fabs(*j) < ACCURACY) ? 0 : (*j)) << '\t';
		}
		if (it + 1 != A.v.end()) out << endl;
	}
	return out;
}

vector<vector<double>> get_vector_from_file(istream& in) {
	vector<vector<double>> v3;
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
	return v3;
}

istream& operator>>(istream& in, Matrix& A) {
	vector<vector<double>> v3;
	string line;
	v3 = get_vector_from_file(in);
	A = Matrix(v3);
	return in;
}

Matrix readbinary(string fname) {
	ifstream in;
	in.open(fname);
	if (!in.is_open()) throw OpenProblem();
	unsigned int M;
	unsigned int N;
	in.read((char* )&M, sizeof(unsigned int));
	in.read((char *)&N, sizeof(unsigned int));
	vector <vector<double>> v3;
	v3.reserve(M);
	for (int i = 0; i < M; i++) {
		vector<double> v1;
		v1.reserve(N);
		for (int j = 0; j < N; j++) {
			double elm;
			in.read((char *)&elm, sizeof(double));
			if (in.eof()) throw FormatProblem();
			v1.push_back(elm);
		}
		v3.push_back(v1);
	}
	in.close();
	return Matrix(v3);
}

void makebinary(const Matrix& A, const string fname = "new_matrix.bin") {
	ofstream out;
	out.open(fname);
	if (!out.is_open()) throw OpenProblem();
	unsigned int M = A.get_rows();
	unsigned int N = A.get_columns();
	out.write((char *)&M, sizeof(unsigned int));
	out.write((char *)&N, sizeof(unsigned int));
	for (auto &col : A.get_matrix()) {
		for (auto &elm : col) {
			out.write((char *)&elm, sizeof(double));
		}
	}
	out.close();
}

Matrix::Matrix(const string& fname) {
	Matrix A;
	if (fname.substr(fname.find_last_of(".") + 1) == "txt") {
		ifstream in;
		in.open(fname);
		if (!in.is_open()) throw OpenProblem();
		in >> A;
		in.close();
	}
	else if (fname.substr(fname.find_last_of(".") + 1) == "bin") {
		A = readbinary(fname);
	}
	else {
		throw UnknownExtention();
	}
	M = A.M;
	N = A.N;
	v.reserve(A.v.size());
	v = A.v;
}


int main() {
	vector<double> v1 = { 0, 1, 2, 5, 7 };
	vector<double> v2 = { 55, 6, 8, 0, 1 };
	vector<double> v3 = { 9, 0, -1, 96, 2 };
	vector<double> v4 = { 0, 4, -3, 343, 3 };
	vector<double> v5 = { 2, 52, 9, -5, 8 };
	vector<vector<double>> v = { v1, v2, v3, v4, v5 };
	vector<double> av1 = { 2, 52, 9, 67, 1, 2, 5, 8, 8, 8, 10 };
	vector<double> av2 = { 2, 52, 9, 67, 1, 40, 5, 16, 16, 16, 20 };
	vector<double> av3 = { 4, 104, 18, 67, 2, 4, 10, 16, 16, 16, 20 };
	vector<double> av4 = { 2, 52, 9, 67, 1, 2, 5, 16, 16, 16, 20 };
	vector<double> av5 = { 2, 52, 9, 67, 1, 41, 5, 16, 16, 16, 20 };
	vector<vector<double>> av = { av1, av2, av3, av4, av5, av1, av3 };

	vector<double> t1 = { 1 };
	vector<double> t2 = { 1, 2 };
	vector<double> t3 = { 1, 2, 3 };
	vector<double> t4 = { 1, 2, 3, 4 };
	vector<double> t5 = { 1, 2, 3, 4, 5 };
	vector<vector<double>> tv3 = { t1, t2, t3, t4 };
	vector<vector<double>> tv2 = { t5, t4, t3, t2, t1 };

	vector<vector<double>> vec01 = { {1}, {2}, {3}, {4}, {5}, {6} };
	vector<vector<double>> vec02 = { {1, 2, 3, 4, 5, 6} };
	vector<vector<double>> vec03 = { { 1, 2, 3, 4, 5 } };
	
	try {
		Matrix A(v);
		IdMatrix B(5);
		UpperTriMatrix C(tv2);
		LowerTriMatrix D(tv3);
		Matrix E(av);
		Matrix vec1(vec01);
		Matrix vec2(vec02);
		Matrix vec3(vec03);
#if LAB == 1:
		cout << A << endl << endl;
		cout << B << endl << endl;
		cout << A + B << endl << endl;
		cout << A - B << endl << endl;
		cout << A * B << endl << endl;
		cout << C << endl << endl;
		cout << D << endl << endl;
		cout << D.Adamar(D) << endl << endl;
		cout << SymmetricMatrix(tv3) << endl << endl;
		cout << "Testing Error" << endl;
		cout << Matrix(tv3) << endl << endl;
		cout << D + A << endl << endl;
#endif
#if LAB == 2:
		cout << A << endl << endl;
		cout << A.Trace() << endl;
		cout << A.Det() << endl;
		cout << scalarMult(vec1, vec2) << endl;
		cout << scalarMult(vec2, vec1) << endl;
		cout << vec1.Norm() << endl;
		cout << vec2.MaxNorm() << endl;
		cout << vec1.MaxNorm() << endl;
		cout << A.Norm() << endl;
		cout << "Testing errors:" << endl;
		//cout << scalarMult(vec1, A) << endl;
		//cout << scalarMult(A, vec2) << endl;
		//cout << scalarMult(A, vec1) << endl;
		//cout << scalarMult(vec1, vec3) << endl;
		//cout << scalarMult(vec3, vec2) << endl;
#endif
#if LAB == 3:
		//cout << AngleVectors(vec2, vec2) << endl << endl;
		//cout << B << endl << endl;
		cout << A << endl << endl;
		//cout << A.Transpose() << endl << endl;
		cout << A.Rank() << endl << endl;
		cout << A.Inverse() << endl << endl;
		cout << A.Inverse() * A << endl << endl;
		cout << A * A.Inverse() << endl << endl;
		//cout << E.Rank() << endl << endl;
		cout << E << endl << endl;
		cout << E.Transpose() << endl << endl;
		cout << "Testing Errors" << endl;
		//cout << AngleVectors(vec1, vec3) << endl << endl;
		//cout << E.Inverse() << endl << endl;
#endif
#if LAB == 4:
		cout << Matrix("txt_actual/scores.txt") << endl << endl;
		cout << A << endl << endl;
		makebinary(A, "sample_binary.bin");
		cout << Matrix("sample_binary.bin") << endl << endl;
		cout << UpperTriMatrix("txt_some/utm.txt") << endl << endl;
		cout << "Testing Errors" << endl;
		cout << Matrix("sample_binary.ababab") << endl << endl;
		cout << Matrix("errorbinary.bin") << endl << endl;
#endif
	}
	catch (MatrixProblems& me) {
		cout << "Matrix ";
		cout << me.what() << endl;
	}
	catch (FileProblem& fe) {
		cout << "File ";
		cout << fe.what() << endl;
	}
	catch (exception& e) {
		cout << e.what() << endl;
	}
	catch (...) {
		cout << "Unknown error" << endl;
	}
}