#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <tuple>
using namespace std;

const double A = -9.3;
const double B = 15;

tuple<double, double, int> divisional(double res, double prevU, double x, int k, double eps) {
	double u = prevU * (x / k);
	if (eps > 1 && k < eps)
		return divisional(u + res, u, x, ++k, eps);
	if (abs(u) > eps)
		return divisional(u + res, u, x, ++k, eps);
	return make_tuple(res, prevU, k);
}

double integer(int x) {
	double res = 1, exp;
	if (x < 0) {
		x = abs(x);
		exp = 1 / M_E;
	}
	else
		exp = M_E;
	for (int i = 0; i < x; i++)
		res *= exp;
	return res;
}

double getX(double left, double right) {
	return (left + right) / 2;
}

int outFirstTable(double _int, double _float) {
	int k;
	cout << "Eps         N delta              R(x)\n";
	for (double eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
		auto divRes = divisional(1, 1, _float, 1, eps);
		if (eps == 1e-8) k = get<2>(divRes);
		printf("%.1e%5i %0.16f %0.16f\n", eps, get<2>(divRes),
			exp(getX(A, B)) - integer(_int) * get<0>(divRes), get<1>(divRes));
	}
	return k;
}

void outSecondTable(int k) {
	double x, _float, _int;
	cout << "X             delta              R(x)\n";
	for (int i = 0; i <= 10; i++) {
		x = A + (B - A) / 10 * i;
		_float = modf(x, &_int);
		auto divRes = divisional(1, 1, _float, 1, k);
		printf("% 8.4f     % 0.16f %0.16f\n", x, abs(exp(x) - integer(_int) * get<0>(divRes)), get<1>(divRes));
	}
}

int main() {
	double num, _int, _float; int n;

	_float = modf(getX(A, B), &_int);
	n = outFirstTable(_int, _float);
	cout << endl;
	outSecondTable(n);

	getchar();
	return 0;
}
