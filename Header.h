#pragma once

#include <iostream>
#include <vector>
#include <iomanip>

#define DEBUG 1
#define DEBUG_LENGTH_STEP_BY_STEP 0
#define DEBUG_DIRECTION_STEP_BY_STEP 0
#define DIMENSION 2

enum methods { STEEPEST_DESCENT, NEWTON };

void printStep(std::vector<double> star, double delta, double a, int k)
{
#if DEBUG
	star[2] = 100 * (star[1] - star[0] * star[0]) + (1 - star[0]) * (1 - star[0]);	// Calculate cost at iteration k+1
	std::cout << std::fixed;
	std::cout << std::setprecision(6);
	std::cout << k << "\t\t" << star[0] << " " << star[1] << "\t\t" << star[2]
		<< "\t\t" << delta << "\t\t" << a << std::endl;
#endif
	return;
}

double L2_norm(std::vector<double> states) 
{
	double norm = 0;
	for (int i = 0; i < DIMENSION; i++)
	{
		norm += states[0] * states[0];
	}
	return sqrt(norm);
}


double eval_f(double x1, double x2)
{
	return 100 * pow( (x2 - (x1 * x1)), 2) + pow((1 - x1), 2);
}

double df_dx1_at_x(double x1, double x2)
{
	return -400 * x1 * (x2 - pow(x1, 2)) - 2 * (1 - x1);
}

double df_dx2_at_x(double x1, double x2)
{
	return 200 * (x2 - x1 * x1);
}

/* Compute Hessian matrix*/
std::vector<std::vector<double>> hessian_f_at_x(std::vector<double> states)
{
	double x1 = states[0];
	double x2 = states[1];

	std::vector<std::vector<double>> H
	 {
		{0, 0},
		{0, 0}
	};

	H[0][0] = -400 * (x2 - x1 * x1) + 800 * x1 * x1 + 2;	// h_11
	H[0][1] = -400 * x1;	// h_12
	H[1][0] = -400 * x1;	// h_21
	H[1][1] = 200;	//h_22

	return H;
}

/* Invert 2 by 2 matrix*/
std::vector<std::vector<double>> invert_2by2(std::vector<std::vector<double>> H)
{
	std::vector<std::vector<double>> H_inverse
	{
		{0, 0},
		{0, 0}
	};

	double det = H[0][0] * H[1][1] - H[0][1] * H[1][0];
	if (det == 0)
	{
		std::cout << "Unable to invert Hessian. Determinant is zero\n";
		return H;
	}
	else
	{
		H_inverse[0][0] = H[1][1] / det;
		H_inverse[0][1] = - H[0][1] / det;
		H_inverse[1][0] = - H[1][0] / det;
		H_inverse[1][1] = H[0][0] / det;

		return H_inverse;
	}
}







