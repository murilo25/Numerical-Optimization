// Murilo Augusto Pinheiro
// 04-Dec-2020
// Solution for Diamond Kinetics - Dynamics Algorithm Implementation.pdf
// main.cpp : Minimizes the Rosenbrock function and displays the result.
// optimization.h : Perform line search optimization to minimize the Rosenbrock function.
// finddirection.h : Calculates the step direction based on method received as an argument. 
//                   (steepest descent and newton's method)
// findsteplength.h : Calculates the step length direction based on the backtracking line search. 
// Header.h : Contains helper functions to evaluate f(x) and its derivatives.


#include "optimization.h"

#define MAX_STEPS 10000
#define TOLERANCE 0.00001
#define INITIAL_CONDITION {1.2,1.2}
#define INITIAL_STEP_LENGTH 1
#define SEARCH_METHOD STEEPEST_DESCENT   //STEEPEST_DESCENT or NEWTON

int main()
{

    std::vector<double> result;
    double start, end;

    start = clock();
    result = optmization(MAX_STEPS, TOLERANCE, INITIAL_CONDITION, INITIAL_STEP_LENGTH, SEARCH_METHOD);
    end = clock();

    std::cout << "\nx1* = " << result[0] << "\nx2* = " << result[1] << "\ncost = " << result[2] << std::endl;

    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    std::cout << "Elapsed time : " << std::fixed
        << time_taken << std::setprecision(5);
    std::cout << " sec " << std::endl;

    return 0;
}
