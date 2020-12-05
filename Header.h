#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>

#define DEBUG 1
#define DIMENSION 2

enum methods { STEEPEST_DESCENT, NEWTON };

std::vector<double> findSearchDirection(std::vector<double> states)
{
	std::vector<double> dir;
	for (int i = 0; i < DIMENSION; i++)
	{
		dir.push_back(0);
		dir.push_back(0);
	}
	//dir[0] = 2 * states[0];	// cost function (x^2 + y^2 + 10) -> derivative: [2x ; 2y]
	//dir[1] = 2 * states[1];

	dir[0] = - 400 * states[0] * (states[1] - pow(states[0],2)) - 2 * (1 - states[0]);
	dir[1] = 200 * (states[1] - states[0] * states[0]);

	return dir;
}

std::vector<double> backtracking()
{
	std::vector<double> step;
	for (int i = 0; i < DIMENSION; i++)
	{
		step.push_back(0);
		step.push_back(0);
	}
	step[0] = 0.002;
	step[1] = 0.002;

	return step;
}


std::vector<double> optmization (int iterMax, double tol, std::vector<double> ic, int method)
{
	std::vector<double> star; // Will store optimization result y*,x*
	star.push_back(ic[0]);	// Store initial condition (x1) in star[0]
	star.push_back(ic[1]);	// Store initial condition (x2) in star[0]
	star.push_back(9999);	// Initialize cost of function as a high value
	double delta = 1;
	double delta_x1 = 1;
	double delta_x2 = 1;
	double k = 0;
	std::vector<double> a; // Initialize step length
	std::vector<double> p; // Initialize search direction

	for (int i = 0; i < DIMENSION; i++) // initialize step and search directions
	{	
		a.push_back(0.1);
		p.push_back(1);
	}

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

				p = findSearchDirection(star);
				a = backtracking();

				for (int i = 0; i < DIMENSION; i++)		// for each dimension
				{
					star[i] = star[i] - a[i] * p[i];	// Steepest descent (gradient descent)
				}

				// Euclidian distance between step k + 1 and k
				delta = abs( sqrt( pow(star[0] - delta_x1,2) + pow(star[1] - delta_x2,2) ));


				k += 1;	// Increment iteration

				if (k > iterMax) 
				{	// Check if it has exceed maximum number of iterations
					delta = -1;
				}

	#if DEBUG
				//star[2] = pow(star[0], 2) + pow(star[1], 2) + 10;	// Calculate cost at iteration k+1
				star[2] = 100 * (star[1] - star[0] * star[0]) + (1 - star[0]) * (1 - star[0]);
				std::cout << std::fixed;
				std::cout << std::setprecision(6);
				std::cout << k << "\t\t\t" << star[0] << " " << star[1] << "\t\t\t" << star[2] 
					<< "\t\t\t" << delta << std::endl;
	#endif

			}
			break;
		case NEWTON:
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
		//star[2] = pow(star[0], 2) + pow(star[1], 2) + 10;	// Calculate cost at iteration k+1 ( cost function (x^2 + y^2 + 10) )
		star[2] = 100 * (star[1] - star[0] * star[0]) + (1 - star[0]) * (1 - star[0]);
	}
	
	return star;
}