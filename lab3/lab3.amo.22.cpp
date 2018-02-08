#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

void outMatrix(vector<vector<double>> mat) {
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j <= mat.size(); j++) {
			cout << mat[i][j] << " ";
		}
		cout << "\n";
	}
}

vector<vector<double>> * methodGaussGordan(vector<vector<double>> * matr) {
	for (int k = 0; k < matr->size(); k++) {
		for (int j = k + 1; j <= matr->size(); j++)
			matr->at(k)[j] = matr->at(k)[j] / matr->at(k)[k];
		matr->at(k)[k] = 1;
		for (int i = 0; i < matr->size(); i++) {
			if (i != k)
				for (int j = k + 1; j <= matr->size(); j++)
					matr->at(i)[j] -= matr->at(k)[j] * matr->at(i)[k];
			if (i != k)
				matr->at(i)[k] = 0;
		}
	}
	return matr;
}

double getGaussSeidelNorm(vector<vector<double>> * matr) {
	vector<double> * maxVect = new vector<double>(matr->size(), 0);

	for (int i = 0; i < matr->size(); i++) {
		for (int j = 0; j < matr->size(); j++) {
			maxVect->at(i) += matr->at(i)[j];
		}
	}
	return * max_element(maxVect->begin(), maxVect->end());
}

vector<double> * methodGaussSeidel(vector<vector<double>> * matr, vector<double> * b, double eps) {
	double norma, q;
	vector<double> * x = new vector<double>(4, 0);
	vector<double> * xTpm = new vector<double>(4, 0);
	vector<double> * normaTmp = new vector<double>(4, 0);

	do {
		norma = 0;
		for (int i = 0; i < matr->size(); i++) {
			x->at(i) = -b->at(i);
			for (int j = 0; j < matr->size(); j++) {
				if (i != j)
					x->at(i) += matr->at(i)[j] * x->at(j);
			}
			x->at(i) = x->at(i) / -matr->at(i)[i];
		}
		for (int i = 0; i < matr->size(); i++) {
			normaTmp->at(i) = abs(x->at(i) - xTpm->at(i));
			xTpm->at(i) = x->at(i);
		}
		norma = * max_element(normaTmp->begin(), normaTmp->end());
		q = getGaussSeidelNorm(matr);
	} while (norma > abs((1 - q) / q) * eps);
	return x;
}
/*
int main() {
	vector<vector<double>> matrGG{
		{ 10, 11, 14, 3, 128 },
		{ 12, 45, 17, 15, 325 },
		{ 9, 16, 31, 5, 240 },
		{ 0, 4, 17, 15, 120 }
	};

	vector<vector<double>> matrGS{
		{ 10, 11, 14, 3 },
		{ 12, 45, 17, 15 },
		{ 9, 16, 31, 5 },
		{ 0, 4, 17, 15 }
	};
	vector<double> b = { 128, 325, 240, 120 };
	vector<double> * x = methodGaussSeidel(&matrGS, &b, 1e-14);

	cout << "Gaus & Gordan method" << endl;
	outMatrix(*(methodGaussGordan(&matrGG)));
	cout << endl << "Gaus & Seidel method" << endl;
	for (int i = 0; i < x->size(); i++) {
		cout << x->at(i) << " ";
	}
	getchar();
	return 0;
}

*/