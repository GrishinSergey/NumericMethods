#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

double myFun(double x) {
	return 10 * x * x * cosh(x) * sin(13 * x);
}

double * getCoefA(int N, double a, double b, double h) {
	double * A = new double[N], x = a;
	for (int i = 0; i < N; x += h, i++)
		A[i] = myFun(x);
	return A;
}

double * getCoefD(double * C, int N, double a, double b, double h) {
	double * D = new double[N];
	D[0] = C[0];
	for (int i = 1; i < N; i++)
		D[i] = (C[i] - C[i - 1]) / h;
	return D;
}

double * getCoefB(double * C, double * D, int N, double a, double b, double h) {
	double * B = new double[N];
	B[0] = 0;
	for (int i = 1; i < N; i++)
		B[i] = (h * C[i] / 2) - (h * h * D[i] / 2) + (myFun(a + h * i) - myFun(a + h * (i - 1))) / h;
	return B;
}

double * sweep(int N, double a, double b, double h) {
	double * C = new double[N],
			* alpha = new double[N],
			* betta = new double[N],
			** slar = new double *[N],
			x1, xi, xn;

	for (int i = 0; i < N; i++)
		slar[i] = new double[N + 1];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N + 1; j++)
			slar[i][j] = 0;

	for (int i = 0; i < N; i++) {
		slar[i][i] = 4 * h;
		if (i > 0)
			slar[i][i - 1] = h;
		if (i < N)
			slar[i][i + 1] = h;
		x1 = a + h * (i - 1);
		xi = a + h * i;
		xn = a + h * (i + 1);
		slar[i][N] = 6 * (myFun(xn) - 2 * myFun(xi) + myFun(x1)) / h;
	}

	alpha[0] = -slar[0][1] / slar[0][0];
	betta[0] =  slar[0][N] / slar[0][0];
	alpha[1] = -slar[1][2] / slar[1][1];
	betta[1] =  slar[1][N] / slar[1][1];

	for (int i = 1; i < N; i++) {
		alpha[i + 1] = -slar[i][i + 1] / (slar[i][i - 1] * alpha[i] + slar[i][i]);
		betta[i + 1] = (slar[i][N] - slar[i][i - 1] * betta[i]) / (slar[i][i - 1] * alpha[i] + slar[i][i]);
	}
	C[N - 1] = (slar[N - 1][N] - slar[N - 1][N - 3] * betta[N - 1]) / (slar[N - 1][N - 3] * alpha[N - 1] + slar[N - 1][N - 2]);
	for (int i = N - 1; i >= 0; i--)
		C[i] = alpha[i + 1] * C[i + 1] + betta[i + 1];
	return C;
}

double getSplainValue(double x, double * A, double * B, double * C, double * D, int N, double a, double b, double h) {
	int i; double xi = a;
	for (i = 0; x >= xi; i++, xi += h);
	xi = a + --i * h;
	return A[i] + B[i] * (x - xi) + C[i] * (x - xi) * (x - xi) / 2 + D[i] * (x - xi) * (x - xi) * (x - xi) / 6;
}

int main() {
	const int N = 50;
	const double A = 0;
	const double B = 1.2;
	const double H = (B - A) / N;
	int count = (int) ((B - A) * N + 1);
	double step = (B - A) / (double)count, x = A + 1e-2, y,
		*a = getCoefA(N, A, B, H),
		*c = sweep(N, A, B, H),
		*d = getCoefD(c, N, A, B, H),
		*b = getCoefB(c, d, N, A, B, H);

	cout << "N = " << N << endl;
	ofstream fout("lab6.csv");

	for (int i = 0; i <= count; i++, x += step) {
		y = getSplainValue(x, a, b, c, d, N, A, B, H);
		fout << x << ";" << y << endl;
	}
	fout.close();
	cout << "finished" << endl;

	system("PAUSE");
	return 0;
}