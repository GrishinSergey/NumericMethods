#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

double f(double x) {	//Main function
	return 2 * sqrt(2 + 3 * sin(x)) - x / 2;
}

double f1(double x) {	//1st derivative
	return (3 * cos(x)) / sqrt(3 * sin(x) + 2) - 0.5;
}

double f2(double x) {	//2nd derivative
	return 2 * (- ((3 * sin(x)) / (2 * sqrt(3 * sin(x) + 2))) - ((9 * cos(x) * cos(x)) / (4 * (3 * sin(x) + 2) * sqrt(3 * sin(x) + 2))));
}

int isolation(double left, double right, double step, double * mas) {
	int i = -1;
	for (double a = left; a <= right; a += step)
		if (f(a) * f(a + step) < 0)
			* (mas + (++i)) = a;
	return i;
}

int * iteration(double * mas, double step, int countRoots) {
	double m, M, a, b, alfa, q, eps, x1, x;
	int n, t = 0, * nIt = (int*) malloc(sizeof(int) * 4);
	for (int i = 0; i <= countRoots; i++) {
		printf(" epr            root%d              accurately\n\n", i + 1);
		a = *(mas + i);
		b = a + step;
		if (abs(f1(a)) > abs(f1(b))) {
			M = f1(a);
			m = f1(b);
		}
		else {
			M = f1(b);
			m = f1(a);
		}
		alfa = 1 / M;
		q = 1 - abs(m / M);
		for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
			x = (b + a) / 2;
			x1 = x;
			n = 0;
			do {
				n++;
				x = x1;
				x1 = x - alfa*f(x);
			} while (abs(x1 - x) > (1 - q) / q * eps);
			if (i == 0) *(nIt + t++) = n;
			printf("%0.e  %20.15f  %20.15f   \n", eps, x1, abs(abs(x1) - abs(x)) * q / (1 - q));
		}
		printf("\n");
	}
	return nIt;
}

int * tangent(double * mas, double step, int countRoots) {
	double a, b, m, eps, x, x1;
	int n, k = 0;
	int * nTan = (int *)malloc(sizeof(int) * 4);
	for (int i = 0; i <= countRoots; i++) {
		printf(" eps            root%d              accurately\n\n", i + 1);
		a = *(mas + i);
		b = a + step;
		m = abs(f1(a));
		for (x1 = a; x1 < b; x1 += 0.001) {
			if (abs(f1(x1)) <= m) m = abs(f1(x1));
		}
		if (f(a) * f2(a) > 0)
			x1 = a;
		else
			x1 = b;
		for (eps = 1e-2; eps >= 1e-14; eps *= 1e-3) {
			x = x1;
			n = 0;
			while (abs(f(x)) / m > eps) {
				n++;
				x -= f(x) / f1(x);
			}
			if (i == 0)
				* (nTan + k++)= n;
			printf("%0.e  %20.15f  %20.15f\n", eps, x, abs(f(x)) / m);
		}
		printf("\n");
	}
	return nTan;
}
/*
int main() {
	double mas[5], step = 0.05, p = 1e+1;
	int * nIt, * nTan, k;

	k = isolation(-10, 10, step, mas);
	nIt = iteration(mas, step, k);
	printf("------------------------------------------------------\n\n");
	nTan = tangent(mas, step, k);
	printf("------------------------------------------------------\n\n");
	printf(" eps   iteration    tangent \n\n");
	for (int i = 0; i <= 4; i++)
		printf("%0.e     %2.0d           %d\n", p /= 1000, nIt[i], nTan[i]);
	system("PAUSE");
}

*/