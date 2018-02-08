#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <functional>

using namespace std;

double simpson(double a, double b, double eps, function<double(double)> fun) {
	int i, n = (int) (1.0 / sqrt(eps)), currentN;
	double h = (b - a) / n, sig1 = 0, sig2 = 0, y0 = fun(a), yn = fun(b),
		currentOdd, xi, currentInt, prevInt, currentEven;
	for (i = 1; i < n; i++)
		if (i % 2 == 0)
			sig1 += fun(a + i * h);
		else
			sig2 += fun(a + i * h);
	prevInt = h / 3 * (4 * sig2 + 2 * sig1 + y0 + yn),
	currentEven = sig2 + sig1;
	currentN = 2 * n;
	do {
		h = (b - a) / currentN;
		currentOdd = 0;
		for (i = 1, xi = a + h; i < currentN; i += 2, xi += 2 * h)
			currentOdd += fun(xi);
		currentInt = h / 3 * (4 * currentOdd + 2 * currentEven + y0 + yn);
		prevInt = currentInt;
		currentEven = currentOdd + currentEven;
		currentN *= 2;
	} while (!(((abs((currentInt - prevInt) / currentInt) / 15) < eps) ||
				((abs(currentInt - prevInt) / 15) < eps)));
	return currentInt;
}

double chebyshev(int n, double x) {
	double tn1, tn = x, tmp = 1;

	if (n == 0) return tmp;
	for (int  i = 1; i < n; i++) {
		tn1 = 2 * x * tn - tmp;
		tmp = tn;
		tn = tn1;
	}
	return tn;
}

double * solve(double ** matrix, int N) {
	double tmp, * result = new double[N];
	for (int j = 0; j < N; j++) {
		tmp = matrix[j][j];
		for (int i = j; i < N + 1; i++)
			matrix[j][i] /= tmp;
		for (int k = j + 1; k < N; k++) {
			tmp = matrix[k][j];
			for (int l = j; l < N + 1; l++)
				matrix[k][l] -= matrix[j][l] * tmp;
		}
	}
	for (int j = N - 1; j >= 0; j--) {
		tmp = 0;
		for (int i = N - 1; i > j; i--)
			tmp += matrix[j][i] * result[i];
		result[j] = matrix[j][N] - tmp;
	}
	return result;
}

double * solveMatrix(double a, double b, double eps, int N, function<double(double)> fun) {
	double ** matrix = new double * [N + 1];
	for (int j = 0; j < N + 1; j++)
		matrix[j] = new double[N + 2];
	for (int j = 0; j < N + 1; j++) {
		for (int i = 0; i < N + 1; i++) {
			if (j != i) matrix[i][j] = simpson(a, b, eps, [&](double x) { return chebyshev(j, x) * chebyshev(i, x); });
						matrix[j][i] = simpson(a, b, eps, [&](double x) { return chebyshev(j, x) * chebyshev(i, x); });
		}
		matrix[j][N + 1] = simpson(a, b, eps, [&](double x) { return chebyshev(j, x) * fun(x); });
	}
	return solve(matrix, N + 1);
}

double chebyshev(double * matrixSolve, int N, double x) {
	double tn1, tn = x, tmp = 1, sum = tmp * matrixSolve[0];

	if (N == 0) return sum;
	sum += tn * matrixSolve[1];
	for (int  i = 1; i < N + 1; i++) {
		tn1 = 2 * x * tn - tmp;
		tmp = tn;
		tn = tn1;
		sum += tn * matrixSolve[i + 1];
	}
	return sum;
}

bool isStandardDeviation(function<double(double)> fun, double * matrixSolve, int N, int a, int b, double eps) {
	return sqrt((simpson(a, b, eps, [&](double x) {
		return (fun(x) - chebyshev(matrixSolve, N, x)) * (fun(x) - chebyshev(matrixSolve, N, x));
	})) / (b - a)) < 0.1;
}

double myFunction(double x) {
	return 10 * x * cosh(x) * sin(13 * x);
}

int main() {

	const double EPS = 1e-7;
	const double X_LEFT = 0, X_RIGHT = 1.2;
	const int ORD_LEFT = 10, ORD_RIGHT = 30;
	ofstream stream;
	double * matrixSolve;

	for (int ord = ORD_LEFT; ord < ORD_RIGHT; ord++) {
		matrixSolve = solveMatrix(X_LEFT, X_RIGHT, EPS, ord, myFunction);
		if (isStandardDeviation(myFunction, matrixSolve, ord, X_LEFT, X_RIGHT, EPS)) {
			stream.open("points.csv");
			for (double x = X_LEFT; x <= X_RIGHT; x += 0.1) {
				cout << setprecision(6) << x << ";" << chebyshev(matrixSolve, ord, x) << endl;
				stream << setprecision(6) << x << ";" << chebyshev(matrixSolve, ord, x) << endl;
			}
			stream.close();
			cout << "N: " << ord << endl;
			break;
		}
	}
	system("PAUSE");
	return 0;

}