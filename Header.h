#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

#define DEBUG 1
#define DIMENSION 2

enum methods { STEEPEST_DESCENT, NEWTON };

double L2_norm(std::vector<double> states) 
{
	return sqrt( (states[0] * states[0]) + (states[1] * states[1]) );
}


double eval_f_at_x(double x1, double x2)
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

std::vector<double> findSearchDirection(std::vector<double> states, int method)
{
	std::vector<double> dir;
	for (int i = 0; i < DIMENSION; i++)
	{
		dir.push_back(0);
		dir.push_back(0);
	}
	std::vector<std::vector<double>> hessianInverse;
	std::vector<double> gradient(2);
	switch (method)
	{
		case STEEPEST_DESCENT:
			/* normalized */
			dir[0] = df_dx1_at_x(states[0], states[1]) / abs(df_dx1_at_x(states[0], states[1]));
			dir[1] = df_dx2_at_x(states[0], states[1]) / abs(df_dx2_at_x(states[0], states[1]));
			/* non normalized */
			//dir[0] = df_dx1_at_x(states[0], states[1]);
			//dir[1] = df_dx2_at_x(states[0], states[1]);
			break;
		case NEWTON:
			// search direction is given by grad_f_at_x * (Hessian_f_at_x)^-1
			gradient[0] = df_dx1_at_x(states[0], states[1]);
			gradient[1] = df_dx2_at_x(states[0], states[1]);
			hessianInverse = invert_2by2(hessian_f_at_x(states));
			dir[0] = hessianInverse[0][0] * gradient[0] + hessianInverse[0][1] * gradient[1];
			dir[1] = hessianInverse[1][0] * gradient[0] + hessianInverse[1][1] * gradient[1];
			break;
		default:
			std::cout << "Unknown method\n";
			break;
	}



	return dir;
}


double findStepLength(std::vector<double> p, std::vector<double> states)
{
	double a;	// step length to be calculated
	double a_bar = 1;
	double rho = 0.8;
	double c = 0.5;
	
	a = a_bar;

	double f_at_x_plus_ap = eval_f_at_x((states[0] + a * p[0]), (states[1] + a * p[1]) );	// f(x_k + a * p_k)
	double f_at_x = eval_f_at_x(states[0], states[1]);	// f(x_k)
	double grad_f_times_pk = df_dx1_at_x(states[0], states[1]) * p[0] + df_dx2_at_x(states[0], states[1]) * p[1];	// transpose(grad_f_x) * p_k
	double secondTerm = c * a * grad_f_times_pk;

	//std::cout << "x1 = " << states[0] << " x2 = " << states[1] << " a = " << a << " p[0] = " << p[0] << " p[1] = " << p[1] << std::endl;
	//std::cout << "x1+a*p0 = " << (states[0] + a * p[0]) << " x2 + a * p1 = " << (states[1] + a * p[1]) << std::endl;
	//std::cout << f_at_x_plus_ap << " <= " << f_at_x << " + " << grad_f_times_pk << " "  << std::endl;
	a = a_bar;
	//while (1);
	
	while ( f_at_x_plus_ap <= f_at_x + secondTerm )	// Backtracking line search condition
	{
		std::cout << states[0] << " " << a << " " << p[0] << " " << states[1] << " " << p[1] << std::endl;

		a = rho * a;
		f_at_x_plus_ap = eval_f_at_x( (states[0] + a * p[0]), (states[1] + a * p[1]) );
		f_at_x = eval_f_at_x(states[0], states[1]);
		grad_f_times_pk = df_dx1_at_x(states[0], states[1]) * p[0] + df_dx2_at_x(states[0], states[1]) * p[1];
		secondTerm = c * a * grad_f_times_pk;
		std::cout << a << "\n";

	}
	

	return a;
}


std::vector<double> optmization (int iterMax, double tol, std::vector<double> ic, int method)
{
	std::vector<double> star; // Will store optimization result x*,f(x*)
	star.push_back(ic[0]);	// Store initial condition (x1) in star[0]
	star.push_back(ic[1]);	// Store initial condition (x2) in star[0]
	star.push_back(9999);	// Initialize cost of function as a high value
	double delta = 1;
	double delta_x1 = 1;
	double delta_x2 = 1;
	double k = 0;	//iterator
	double a = 1; // Initialize step length
	std::vector<double> p; // declare search direction vector

#if DEBUG
	std::cout << "Iter: \t\t\t\t x1 & x2: \t\t\t\t cost: \t\t\t\t delta: \n";
#endif
	switch (method) 
	{
		case STEEPEST_DESCENT:	// create function to improve code readability
			while (delta > tol) 
			{
				delta_x1 = star[0];
				delta_x2 = star[1];
				
				p = findSearchDirection(star,STEEPEST_DESCENT);
				a = findStepLength(p, star);

				for (int i = 0; i < DIMENSION; i++)		// for each dimension
				{
					star[i] = star[i] - a * p[i];	// Steepest descent (gradient descent)
				}

				// Euclidian distance between step k + 1 and k
				delta = abs( sqrt( pow(star[0] - delta_x1,2) + pow(star[1] - delta_x2,2) ));


				k += 1;	// Increment iteration

				if (k > iterMax) 
				{	// Check if it has exceed maximum number of iterations
					delta = -1;
				}

	#if DEBUG
				star[2] = 100 * (star[1] - star[0] * star[0]) + (1 - star[0]) * (1 - star[0]);	// Calculate cost at iteration k+1
				std::cout << std::fixed;
				std::cout << std::setprecision(6);
				std::cout << k << "\t\t\t" << star[0] << " " << star[1] << "\t\t\t" << star[2] 
					<< "\t\t\t" << delta << std::endl;
	#endif

			}
			break;
		case NEWTON:
			while(delta > tol)
			{ 
				double delta_x = star[0];
				double delta_y = star[1];

				p = findSearchDirection(star,NEWTON);
				a = findStepLength(p, star);

				for (int i = 0; i < DIMENSION; i++)		// for each dimension
				{
					star[i] = star[i] - a * p[i];	// Steepest descent (gradient descent)
				}

				// Euclidian distance between step k + 1 and k
				delta = abs(sqrt(pow(star[0] - delta_x1, 2) + pow(star[1] - delta_x2, 2)));

				k += 1;

				if (k > iterMax)
				{	// Check if it has exceed maximum number of iterations
					delta = -1;
				}

#if DEBUG
				star[2] = 100 * (star[1] - star[0] * star[0]) + (1 - star[0]) * (1 - star[0]);	// Calculate cost at iteration k+1
				std::cout << std::fixed;
				std::cout << std::setprecision(6);
				std::cout << k << "\t\t\t" << star[0] << " " << star[1] << "\t\t\t" << star[2]
					<< "\t\t\t" << delta << std::endl;
#endif
			}
			break;
		default:
			std::cout << "Unknown method\n";
			break;
		}

	if (delta > tol || k > iterMax)		// if optimization did not converge, set result to -1 
	{	
		star[0] = -1;
		star[1] = -1;
		star[2] = -1;
		std::cout << "Optimization did not converge\n";
	}
	else	// if converged, calculate cost
	{	
		std::cout << "Optmization converged\n";
		star[2] = 100 * (star[1] - star[0] * star[0]) + (1 - star[0]) * (1 - star[0]);	// calculate f(x*)
	}
	
	return star;
}