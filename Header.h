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

/* Find search direction */
std::vector<double> findSearchDirection(std::vector<double> states, int method)
{
	std::vector<double> dir(2);
	std::vector<double> gradient(2);
	std::vector<std::vector<double>> hessianInverse;
	switch (method)
	{
		case STEEPEST_DESCENT:
			/* normalized */
			//dir[0] = df_dx1_at_x(states[0], states[1]) / abs(df_dx1_at_x(states[0], states[1]));
			//dir[1] = df_dx2_at_x(states[0], states[1]) / abs(df_dx2_at_x(states[0], states[1]));
			/* non normalized */
			dir[0] = df_dx1_at_x(states[0], states[1]);
			dir[1] = df_dx2_at_x(states[0], states[1]);
			break;
		case NEWTON:
			// search direction is given by (Hessian_f_at_x)^-1 * grad_f_at_x
			gradient[0] = df_dx1_at_x(states[0], states[1]);
			gradient[1] = df_dx2_at_x(states[0], states[1]);
			hessianInverse = invert_2by2(hessian_f_at_x(states));
			// Compute (Hessian_f_at_x)^-1 * grad_f_at_x
			dir[0] = hessianInverse[0][0] * gradient[0] + hessianInverse[0][1] * gradient[1];
			dir[1] = hessianInverse[1][0] * gradient[0] + hessianInverse[1][1] * gradient[1];
			/* normalize */
			//dir[0] = dir[0] / abs(dir[0]);
			//dir[1] = dir[1] / abs(dir[1]);
			break;
		default:
			std::cout << "Unknown method\n";
			break;
	}
	//multiply by -1 to convert gradient from ascent to descent
	dir[0] = -dir[0];
	dir[1] = -dir[1];

	#if DEBUG_DIRECTION_STEP_BY_STEP
		std::cout << dir[0] << "\t" << dir[1] << std::endl;
	#endif
	return dir;
}

/* Find step length using backtracking line search */
double findStepLength(std::vector<double> p, std::vector<double> states, double a_bar)
{	
	double a;	// step length to be calculated
	double rho = 0.5;
	double c = 0.15;
	
	a = a_bar;

	// Calculate terms for backtracking search
	double f_at_x_plus_ap = eval_f((states[0] + a * p[0]), (states[1] + a * p[1]) );	// f(x_k + a * p_k)
	double f_at_x = eval_f(states[0], states[1]);	// f(x_k)
	double df_dx1_xk = df_dx1_at_x(states[0], states[1]);	// df_dx1(x_k)
	double df_dx2_xk = df_dx2_at_x(states[0], states[1]);	// df_dx2(x_k)
	double grad_f_times_pk = df_dx1_xk * p[0] + df_dx2_xk * p[1];	// transpose(grad_f_x) * p_k
	double secondTerm = c * a * grad_f_times_pk;
	
	#if DEBUG_LENGTH_STEP_BY_STEP
		std::cout << "F(x+ap)\t\t" << "F(x)\t\t" << "secondTerm\t\t" << "step" << std::endl;
	#endif

	while ( f_at_x_plus_ap > f_at_x + secondTerm )	// Backtracking line search condition
	{
		#if DEBUG_LENGTH_STEP_BY_STEP
			std::cout << f_at_x_plus_ap << "\t\t" << f_at_x << "\t\t" << secondTerm << "\t\t" << a << std::endl;
		#endif
		a = rho * a;
		// Update terms for convergence check
		f_at_x_plus_ap = eval_f( (states[0] + a * p[0]), (states[1] + a * p[1]) );
		secondTerm = c * a * grad_f_times_pk;

	}
	
	return a;
}


std::vector<double> optmization (int iterMax, double tol, std::vector<double> ic, double a0, int method)
{
	std::vector<double> star; // Declare variable to store optimization result x1*,x2*,f(x*)
	for (int i = 0; i < DIMENSION; i++)	// Initialize star variable
	{
		star.push_back(ic[i]);
	}
	star.push_back(eval_f(ic[0], ic[1]));

	double delta = 1;
	double previous_x1 = 1;
	double previous_x2 = 1;
	int k = 0;	// iterator
	double a;	// step length

	std::vector<double> p; // Declare search direction vector
	for (int i = 0; i < DIMENSION; i++)	// Initialize p variable
	{
		p.push_back(0);
	}
	
	#if DEBUG // Display header
		std::cout << "Iter: \t\t\t x1 & x2: \t\t cost: \t\t\t delta: \t\t step: \n";
	#endif
	#if DEBUG_DIRECTION_STEP_BY_STEP // Display header
		std::cout << "dir[0]:  dir[1]: \n";
	#endif

	while (delta > tol)		// Optimization loop
	{
		//std::cout << "iter " << k << "\n";
		previous_x1 = star[0];	// store current x1 to compute delta after one step
		previous_x2 = star[1]; // store current x2 to compute delta after one step

		p = findSearchDirection(star, method);
		a = findStepLength(p, star, a0);

		for (int i = 0; i < DIMENSION; i++)		
		{
			star[i] = star[i] + a * p[i];	// take a step of length a and direction p[i] for each dimension
		}

		// Euclidian distance between step k + 1 and k to check if solution converged
		delta = abs(sqrt(pow(star[0] - previous_x1, 2) + pow(star[1] - previous_x2, 2)));

		k += 1;	// Increment iteration

		if (k > iterMax)	// Check if it has exceed maximum number of iterations
		{	
			delta = -1;
		}

		printStep(star, delta, a, k);	// Display results after 1 iteration

	}

	if (delta > tol || k > iterMax)		// If optimization did not converge
	{	
		star[0] = -1;
		star[1] = -1;
		star[2] = -1;
		std::cout << "Optimization did not converge\n";
	}
	else	// If converged, calculate cost
	{	
		std::cout << "Optmization converged\n";
		star[2] = eval_f(star[0],star[1]);	// Calculate f(x*)
	}
	
	return star;
}

