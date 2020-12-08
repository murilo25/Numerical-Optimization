#pragma once

#include "Header.h"

/* Find step length using backtracking line search */
double findStepLength(std::vector<double> p, std::vector<double> states, double a_bar)
{
	double a;	// step length to be calculated
	double rho = 0.5;
	double c = 0.15;

	a = a_bar;

	// Calculate terms for backtracking search
	double f_at_x_plus_ap = eval_f((states[0] + a * p[0]), (states[1] + a * p[1]));	// f(x_k + a * p_k)
	double f_at_x = eval_f(states[0], states[1]);	// f(x_k)
	double df_dx1_xk = df_dx1_at_x(states[0], states[1]);	// df_dx1(x_k)
	double df_dx2_xk = df_dx2_at_x(states[0], states[1]);	// df_dx2(x_k)
	double grad_f_times_pk = df_dx1_xk * p[0] + df_dx2_xk * p[1];	// transpose(grad_f_x) * p_k
	double secondTerm = c * a * grad_f_times_pk;

#if DEBUG_LENGTH_STEP_BY_STEP
	std::cout << "F(x+ap)\t\t" << "F(x)\t\t" << "secondTerm\t\t" << "step" << std::endl;
#endif

	while (f_at_x_plus_ap > f_at_x + secondTerm)	// Backtracking line search condition
	{
#if DEBUG_LENGTH_STEP_BY_STEP
		std::cout << f_at_x_plus_ap << "\t\t" << f_at_x << "\t\t" << secondTerm << "\t\t" << a << std::endl;
#endif
		a = rho * a;
		// Update terms for convergence check
		f_at_x_plus_ap = eval_f((states[0] + a * p[0]), (states[1] + a * p[1]));
		secondTerm = c * a * grad_f_times_pk;

	}

	return a;
}