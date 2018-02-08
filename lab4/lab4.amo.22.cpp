#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
using namespace std;

double integral(double x) {
	return cos(x) * (1 + 2 * pow(M_E, sin(x)));
}

double antiderative(double x) {
	return 2 * pow(M_E, sin(x)) + sin(x);
}

double trapezium(double a, double b, double h, double (*f) (double)) {
	double sum = (f(a) + f(b)) / 2, xi;
	int n = (int) ceil((b - a) / h);
	h = (b - a) / n;
	xi = a + h;
	for (int i = 1; i < n; i++) {
		sum += f(xi);
		xi += h;
	}
	return h * sum;
}

tuple<long double, long double> runge(double a, double b, double eps, double (*f) (double)) {
	int n = (int) ceil(1 / sqrt(eps));
	double val1, val2;
	do {
		val1 = trapezium(a, b, (b - a) / n, f);
		val2 = trapezium(a, b, (b - a) / (2 * n), f);
		n *= 2;
	} while (abs(val1 - val2) / 3 > eps);
	return make_tuple(val2, (b - a) / n);
}

int main() {
	const double a = 0.0;
	const double b = 15.0;
	const vector<double> steps = { 0.0957658, 0.00957658, 0.000957658, 0.0000957658 };
	vector<double> errors;
	double value, exactValue, error, eps = 0.1;
	int i = 0;

	exactValue = antiderative(b) - antiderative(a);

	printf("trapezium:\n eps\t\t step\t\t exact value\t     error\n");
	for (double eps = 0.1; eps > 1e-7; eps *= 1e-2) {
		value = trapezium(a, b, steps[i], &integral);
		error = abs(value - exactValue);
		printf("%.12f  %.12f  %.15f  %.16f\n", eps, steps[i], exactValue, error);
		errors.push_back(error);
		i++;
	}

	printf("\n\nrunge:\n expected error\t     step\t     prodused error\n");
	for (int i = 0; i < (int) errors.size(); i++) {
		auto res = runge(a, b, errors[i], integral);
		printf("%.16f  %.12f  %.16f\n", errors[i], get<1>(res), abs(get<0>(res) - exactValue));
		eps *= 0.01;
	}
	system("PAUSE");
	return 0;
}