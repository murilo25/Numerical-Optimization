#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

#define DEBUG 0

enum methods { STEEPEST_DESCENT, NEWTON };

std::vector<double> optmization (int iterMax, double tol, std::vector<double> ic, int method)
{
	std::vector<double> star; // Will store optimization result y*,x*
	star.push_back(ic[0]);	// Store initial condition (x1) in star[0]
	star.push_back(ic[1]);	// Store initial condition (x2) in star[0]
	star.push_back(9999);	// Initialize cost of function as a high value
	double delta = 1;
	double k = 0;
	double a = 0.1; // Step length
	double p; // Search direction

	switch (method) {
	case STEEPEST_DESCENT:	// create function to improve code readability
#if DEBUG
			std::cout << "Iter: \t\t\t\t x1: \t\t\t\t cost: \t\t\t\t delta: \n";
#endif
		while (delta > tol) {

			delta = star[0];	// store x_k to calculate delta after 1 iteration
			p = 2 * star[0] + 7;	// cost function (x^2 + 7x +10) -> derivative: 2x + 7
			star[0] = star[0] - a * p;	// Steepest descent (gradient descent)
			delta = abs((star[0] - delta) / delta);	// Check for convergence

			k += 1;	// Increment iteration

			if (k > iterMax) {	// Check if it has exceed maximum number of iterations
				delta = -1;
			}

#if DEBUG
			star[2] = pow(star[0], 2) + 7 * star[0] + 10;	// Calculate cost at iteration k+1
			std::cout << std::fixed;
			std::cout << std::setprecision(6);
			std::cout << k << "\t\t\t" << star[0] << "\t\t\t" << star[2] << "\t\t\t" << delta << std::endl;
#endif

		}
		break;
	case NEWTON:
		break;
	default:
		std::cout << "Unknown method\n";
		break;
	}

	if (delta > tol || k > iterMax) {	// if optimization did not converge, set result to -1 
		star[0] = -1;
		star[1] = -1;
		star[2] = -1;
		std::cout << "Optimization did not converge\n";
	}
	else {	// if converged, calculate cost
		std::cout << "Optmization converged\n";
		star[2] = pow(star[0], 2) + 7 * star[0] + 10;
	}
	
	return star;
}