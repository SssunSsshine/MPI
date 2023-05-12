#include  <mpi.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

void input(double& x1, double& xn, double& eps, int& n, int procnum)
{
	std::string str;
	do {
		std::cout << "Input first value: ";
		std::cin >> str;
		std::istringstream istr(str);
		if (!(istr >> x1)) {
			x1 = -3;
		}
	} while (x1 >= 1 || x1 <= -1);

	std::cout << "\n";
	do {
		std::cout << "Input last value: ";
		std::cin >> str;
		std::istringstream istr(str);
		if (!(istr >> xn)) {
			xn = -3;
		}
	} while (xn >= 1 || xn <= -1 || xn < x1);
	std::cout << "\n";

	do {
		std::cout << "Input epsilon: ";
		std::cin >> str;
		std::istringstream istr(str);
		if (!(istr >> eps)) {
			eps = -3;
		}
	} while (eps >= 1 || eps <= 0);
	std::cout << "\n";

	do {
		std::cout << "Input n: ";
		std::cin >> str;
		std::istringstream istr(str);
		if (!(istr >> n)) {
			n = -3;
		}

	} while (n <= 0 || n < procnum);
	std::cout << "\n";
}

void print_x(int n, double* globaldata)
{
	std::cout
		<< std::left
		<< std::setw(20)
		<< "|X Values";
	for (int i = 0; i < n; i++)
		std::cout
		<< std::left
		<< std::setw(10)
		<< "|" + std::to_string(globaldata[i]);
	std::cout << std::endl;
}

void print_eps_result(int n, double* resglobal)
{
	std::cout
		<< std::left
		<< std::setw(20)
		<< "|epsilon results";
	for (int i = 0; i < n; i++)
		std::cout
		<< std::left
		<< std::setw(10)
		<< "|" + std::to_string(resglobal[i]);
	std::cout << std::endl;
}

void print_result_values(int n, double* res)
{
	std::cout
		<< std::left
		<< std::setw(20)
		<< "|Results";
	for (int i = 0; i < n; i++)
		std::cout
		<< std::left
		<< std::setw(10)
		<< "|" + std::to_string(res[i]);
	std::cout << std::endl;
}

void print_header(int n)
{
	std::cout
		<< std::left
		<< std::setw(20)
		<< "|";
	for (int i = 1; i <= n; i++) {
		std::cout
			<< std::left
			<< std::setw(10)
			<< "|X" + std::to_string(i);
	}
	std::cout << std::endl;
}

void bord(int n)
{
	for (int i = 0; i < (10 + n) + 10 * n; i++)
	{
		std::cout << "_";
	}
	std::cout << std::endl;
}

void print_table(int n, double eps, double* globaldata, double* resglobal, double* res)
{
	bord(n);
	print_header(n);
	bord(n);

	print_x(n, globaldata);
	bord(n);

	print_eps_result(n, resglobal);
	bord(n);

	print_result_values(n, res);
	bord(n);
}

//вычисление ряда тэйлора для числа x с точностью eps
double calculate(double x, double eps)
{
	double cur = x/2;
	double res = 1 + 0.5 * x;
	int i = 1;
	while (abs(cur) >= eps) {
		double num = x * (2 * i - 1);
		double de = 2 * (i + 1);
		cur *= - num / de;
		res += cur;
		i++;
	}
	return res;
}

int main(int argc, char* argv[])
{
	int procnum;
	int myrank;
	double x1;
	double xn;
	int n;
	double eps;
	double* globaldata = NULL;
	double* localdata = NULL;
	double* reslocal = NULL;
	double* resglobal = NULL;
	double* res = NULL;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &procnum);

	if (myrank == 0) {
		input(x1, xn, eps, n, procnum);

		globaldata = new double[n];
		res = new double[n];
		resglobal = new double[n];

		double step = (xn - x1) / n;
		double tmp = x1;
		for (int i = 0; i < n; i++) {
			globaldata[i] = tmp;
			res[i] = sqrt(1 + tmp);
			tmp = tmp + step;
		}
	}

	MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	double k = ceil((double)n / (double)procnum);
	localdata = new double[k];
	reslocal = new double[k];

	MPI_Scatter(globaldata, k, MPI_DOUBLE, localdata, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = 0; i < k; i++) {
		reslocal[i] = calculate(localdata[i], eps);
	}

	MPI_Gather(reslocal, k, MPI_DOUBLE, resglobal, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (myrank == 0) {
		print_table(n, eps, globaldata, resglobal, res);
	}

	MPI_Finalize();
	return 0;
}

