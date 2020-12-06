// Murilo Augusto Pinheiro
// 04-Dec-2020
// Solution for Diamond Kinetics - Dynamics Algorithm Implementation.pdf
// SteepestDescent.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Header.h : This files contains all functions and includes used in this file.



#include "Header.h"

#define MAX_STEPS 9000
#define TOLERANCE 0.00001
#define INITIAL_CONDITION {1.2,1.2}
#define INITIAL_STEP_LENGTH 1
#define SEARCH_METHOD STEEPEST_DESCENT   //STEEPEST_DESCENT or NEWTON

int main()
{
    std::vector<double> result;

    result = optmization(MAX_STEPS, TOLERANCE, INITIAL_CONDITION, INITIAL_STEP_LENGTH, SEARCH_METHOD);

    std::cout << "\nx1* = " << result[0] << "\nx2* = " << result[1] << "\ncost = " << result[2] << std::endl;

    return 0;
}
