#pragma once

#include "findstep.h"
#include "finddirection.h"

/* Minimize the Rosenbrock function using a line search approach */
std::vector<double> optmization(int iterMax, double tol, std::vector<double> ic, double a0, int method)
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
		star[2] = eval_f(star[0], star[1]);	// Calculate f(x*)
	}

	return star;
}